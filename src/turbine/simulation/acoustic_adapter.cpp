#include "turbine/acoustic_adapter.hpp"
#include "lookup_diagnostics.hpp"
namespace turbine {
AcousticConfiguration::AcousticConfiguration(const Case &c) {
    parameters = aeroacoustics::read_aa_input(c.acoustic.path.string()).first;
    parameters.airdens = c.rho;
    parameters.kinvisc = c.nu;
    parameters.spdsound = c.sound_speed;
    observers = aeroacoustics::read_observers(c.acoustic.file("ObserverLocations").string());

    for (const auto &s : c.stations)
        span.push_back(s.span);

    aeroacoustics::validate(parameters);
    sample_interval = c.acoustic.number("DT_AA", c.dt);
    start = c.acoustic.number("AAStart");
    blade_percent = c.acoustic.number("BldPrcnt");
    hub_height = c.structure.number("TowerHt") + c.structure.number("Twr2Shft");
    ti_method = c.acoustic.integer("TICalcMeth");
}
aeroacoustics::AcousticDriver AcousticConfiguration::make_driver() const {
    return {parameters, span,          blades,     observers, sample_interval,
            start,      blade_percent, hub_height, ti_method};
}
AcousticInputAdapter::AcousticInputAdapter(const Case &c, const AcousticConfiguration &config,
                                           std::shared_ptr<const SurfaceSet> surfaces)
    : nodes_(config.blades), tables_(c.stations.size()), viscosity_(c.nu), surfaces_(std::move(surfaces)) {
    if (config.blades != 3)
        throw std::invalid_argument("Acoustic adapter requires three blades");
    const auto &parameters = config.parameters;
    for (const auto &station : c.stations)
        airfoil_ids_.push_back(station.airfoil);
    for (int b = 0; b < 3; ++b)
        for (std::size_t j = 0; j < c.stations.size(); ++j) {
            const auto &s = c.stations[j];
            const auto &af = c.airfoils[s.airfoil];
            const auto *surface = surfaces_ ? surfaces_->find(s.airfoil) : nullptr;
            aeroacoustics::Node n;
            n.section.chord = s.chord;
            n.section.stall_deg = af.input.number("alpha1");
            n.airfoil_reference = af.reference;
            if (parameters.timod == 2) {
                const auto thickness = aeroacoustics::guidati_thickness(af.coordinates);
                n.section.thickness_1p = thickness[0];
                n.section.thickness_10p = thickness[1];
            }
            if (surface) {
                n.section.tabulated_boundary_layer = true;
                n.section.te_angle = surface->te_angle_deg;
                n.section.te_thickness = surface->te_thickness_m;
            } else if (parameters.bluntmod) {
                turbine::InputFile bl(af.input.file("BL_file"));
                n.section.te_angle = bl.number("TEAngle");
                n.section.te_thickness = bl.number("TEThick");
            }
            if (!surface && (parameters.x_blmethod == 2 || parameters.tbltemod == 2) && b == 0)
                tables_[j].emplace(aeroacoustics::BLTable::read(af.input.file("BL_file").string()));
            nodes_[b].push_back(n);
        }
}
const std::vector<std::vector<aeroacoustics::Node>> &
AcousticInputAdapter::update(const RotorOutput &aerodynamic, double time, bool sampling,
                             std::size_t first_node) {
    if (nodes_.empty() || first_node >= nodes_.front().size())
        throw std::invalid_argument("Invalid acoustic node selection");
    for (std::size_t b = 0; b < nodes_.size(); ++b)
        if (aerodynamic.blades[b].size() != nodes_[b].size())
            throw std::invalid_argument("Acoustic adapter node count mismatch");
    for (int b = 0; b < 3; ++b)
        for (std::size_t j = 0; j < nodes_[b].size(); ++j) {
            const auto &a = aerodynamic.blades[b][j];
            auto &node = nodes_[b][j];
            node.aero_center = a.motion.position;
            node.inflow = a.wind;
            node.section.speed = a.speed;
            node.section.alpha_deg = a.alpha / turbine::deg;
            for (int i = 0; i < 3; ++i)
                for (int k = 0; k < 3; ++k)
                    node.global_to_local[3 * i + k] = a.motion.orientation[i][k];
            const auto *surface = surfaces_ ? surfaces_->find(airfoil_ids_[j]) : nullptr;
            if (surface && (surface->polar || (sampling && j >= first_node))) {
                // Validate active polar states even on non-acoustic steps. BEM trial angles
                // retain the complete supplied polar; coverage concerns physical states.
                const auto bl = surface->at(node.section.alpha_deg, a.speed * node.section.chord / viscosity_,
                                            node.section.chord);
                if (sampling && j >= first_node)
                    node.section.bl = bl;
            } else if (sampling && j >= first_node && tables_[j]) {
                diagnostics::LookupLocation location({time, std::size_t(b + 1), j + 1, "BL_interpolate"});
                node.section.bl = tables_[j]->interpolate(
                    node.section.alpha_deg, a.speed * node.section.chord / viscosity_, node.section.chord);
            }
        }

    return nodes_;
}
} // namespace turbine
