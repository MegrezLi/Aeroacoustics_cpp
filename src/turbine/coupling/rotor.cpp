#include "turbine/rotor.hpp"
#include "lookup_diagnostics.hpp"
namespace turbine {
Rotor::Rotor(const Case &c) : Rotor(TurbineModel(c)) {}
Rotor::Rotor(TurbineModel model) : model_(std::move(model)), structure_(model_.data()) {
    const auto &c = model_.data();
    options_.tip_loss = c.aero.flag("TipLoss");
    options_.hub_loss = c.aero.flag("HubLoss");
    options_.tangential = c.aero.flag("TanInd");
    options_.axial_drag = c.aero.flag("AIDrag");
    options_.tangential_drag = c.aero.flag("TIDrag");
    options_.tolerance = c.aero.number("IndToler", 5e-10);
    options_.max_iterations = c.aero.integer("MaxIter");
    skew_ = {c.aero.integer("Skew_Mod") == 1 && c.aero.number("SkewRedistr_Mod", 1) == 1,
             c.aero.number("SkewRedistrFactor", 15 * pi / 32)};
    BladeStructure reference = structure_;
    reference.pitch = {0, 0, 0};
    reference.initial_azimuth -= c.structure.number("Azimuth") * deg;
    std::vector<double> arc;
    Vec3 last{};
    double distance = structure_.hub_radius;
    for (const auto &s : c.stations) {
        const Vec3 p{s.curve, s.sweep, s.span};
        distance += norm(p - last);
        arc.push_back(distance);
        last = p;
    }
    for (double z : arc) {
        tip_constant_.push_back(1.5 * (arc.back() - z) / z);
        hub_constant_.push_back(1.5 * (z - structure_.hub_radius) / structure_.hub_radius);
    }
    for (int b = 0; b < 3; ++b) {
        std::vector<ReferenceNode> sr, ar;
        std::vector<Vec3> sp, ap;
        for (std::size_t j = 0; j < reference.nodes.size(); ++j) {
            const auto m = reference.motion(0, b, {}, reference.nodes[j]);
            sr.push_back({m.position, m.orientation});
            if (j > 0 && j + 1 < reference.nodes.size())
                sp.push_back(m.position);
        }
        const auto root = reference.blade_basis(0, b);
        const Vec3 origin = reference.hub + reference.hub_radius * root[2];
        for (const auto &s : c.stations) {
            ar.push_back({origin + multiply(transpose(root), Vec3{s.curve, s.sweep, s.span}),
                          multiply(euler_matrix({0, s.curve_angle, -s.twist}), root)});
            ap.push_back(ar.back().position);
            previous_phi_[b].push_back(0);
            airfoils_[b].push_back(UnsteadyAirfoil(model_.airfoil(s.airfoil), s.chord, c.dt, c.sound_speed));
        }
        motion_maps_.emplace_back(sr, ar);
        load_maps_.emplace_back(ap, sp);
    }
}
RotorOutput Rotor::evaluate(double time, const RotorState &state) const {
    RotorOutput y;
    RotorWorkspace workspace;
    evaluate_into(time, state, y, workspace);
    return y;
}
void Rotor::evaluate_into(double time, const RotorState &state, RotorOutput &y,
                          RotorWorkspace &workspace) const {
    const auto &c = model_.data();
    y.average_velocity = {};
    std::array<Matrix3, 3> basis;
    for (int b = 0; b < 3; ++b)
        basis[b] = structure_.blade_basis(time, b);
    const std::size_t count = c.stations.size();
    for (int b = 0; b < 3; ++b) {
        structure_.motions_into(basis[b], state[b], workspace.structural);
        motion_maps_[b].transfer_into(workspace.structural, workspace.aerodynamic);
        const auto &aerodynamic = workspace.aerodynamic;
        y.blades[b].resize(count);
        // Removing prescribed pitch from the blade basis is equivalent to the
        // Euler decomposition in Calculate_MeshOrientation_NoSweepPitchTwist.
        const auto &root = basis[b];
        const double pitch = structure_.pitch[b];
        const Matrix3 unpitched{{std::cos(pitch) * root[0] + std::sin(pitch) * root[1],
                                 -std::sin(pitch) * root[0] + std::cos(pitch) * root[1], root[2]}};
        for (std::size_t j = 0; j < count; ++j) {
            auto &a = y.blades[b][j];
            const auto &s = c.stations[j];
            a.motion = aerodynamic[j];
            a.wind = c.wind.at(a.motion.position);
            y.average_velocity = y.average_velocity + a.wind - a.motion.velocity;
            const auto angles = euler_angles(multiply(a.motion.orientation, transpose(unpitched)));
            a.annulus = multiply(euler_matrix({0, angles[1], 0}), unpitched);
            const Vec3 relative = a.wind - a.motion.velocity,
                       hub_relative = a.motion.position - structure_.hub;
            a.bem = {norm(hub_relative - dot(hub_relative, structure_.shaft) * structure_.shaft),
                     s.chord,
                     -angles[2],
                     dot(relative, a.annulus[0]),
                     dot(relative, a.annulus[1]),
                     hub_constant_[j],
                     tip_constant_[j]};
        }
    }
    y.average_velocity = y.average_velocity / (3 * count);
    const double axial_velocity = dot(y.average_velocity, structure_.shaft);
    y.skew = norm(y.average_velocity) < 1e-14
                 ? 0
                 : std::acos(std::clamp(axial_velocity / norm(y.average_velocity), -1., 1.));
    const Vec3 transverse = axial_velocity * structure_.shaft - y.average_velocity;
    Vec3 disk_y{}, disk_z{};
    if (norm(transverse) > 1e-14) {
        disk_y = unit(transverse);
        disk_z = cross(y.average_velocity, structure_.shaft) / norm(transverse);
    } else {
        const auto &base = basis[0];
        disk_y = unit(cross(base[2], structure_.shaft));
        disk_z = cross(structure_.shaft, disk_y);
    }
    double max_radius = 0;
    for (const auto &blade : y.blades)
        for (const auto &a : blade)
            max_radius = std::max(max_radius, a.bem.radius);
    const double tsr =
        std::abs(axial_velocity) < 1e-14 ? 1e10 : std::abs(structure_.omega * max_radius / axial_velocity);
    const double weight = tsr >= 2 ? 1 : tsr <= 1 ? 0 : .5 * (1 - std::cos(pi * (tsr - 1)));
    for (int b = 0; b < 3; ++b) {
        const auto &root = basis[b];
        const double azimuth = std::atan2(-dot(root[2], disk_y), dot(root[2], disk_z));
        double radius = 0;
        for (const auto &a : y.blades[b])
            radius = std::max(radius, a.bem.radius);
        for (std::size_t j = 0; j < count; ++j) {
            auto &a = y.blades[b][j];
            const auto &s = c.stations[j];
            const auto induction = [&] {
                diagnostics::LookupLocation location({time, std::size_t(b + 1), j + 1, "BEM_trial"});
                return solve_bem(options_, a.bem, c.airfoils[s.airfoil], previous_phi_[b][j]);
            }();
            a.root_phi = induction.phi;
            a.axial = std::clamp(induction.axial, -1., 1.5);
            a.tangential = std::clamp(induction.tangential, -1., 1.);
            const bool fixed =
                (options_.tip_loss && tip_constant_[j] == 0) || (options_.hub_loss && hub_constant_[j] == 0);
            if (!fixed && skew_.redistribute)
                a.axial = skew_axial(a.axial, y.skew, a.bem.radius / radius, azimuth, skew_.factor);
            a.axial *= weight;
            a.tangential *= weight;
            const double vx = a.bem.vx * (1 - a.axial), vy = a.bem.vy * (1 + a.tangential);
            a.phi = std::atan2(vx, vy);
            a.alpha = std::remainder(a.phi - a.bem.twist, 2 * pi);
            a.speed = std::hypot(vx, vy);
            {
                diagnostics::LookupLocation location({time, std::size_t(b + 1), j + 1, "UA_evaluate"});
                a.coefficients = airfoils_[b][j].evaluate(a.alpha, a.speed);
            }
            const auto cf = a.coefficients;
            const double q = .5 * c.rho * a.speed * a.speed;
            a.load.force =
                multiply(transpose(a.annulus),
                         Vec3{(cf.cl * std::cos(a.phi) + cf.cd * std::sin(a.phi)) * q * s.chord,
                              -(cf.cl * std::sin(a.phi) - cf.cd * std::cos(a.phi)) * q * s.chord, 0});
            a.load.moment = (cf.cm * q * s.chord * s.chord) * a.annulus[2];
        }
    }
}
void Rotor::advance_airfoils(const RotorOutput &y, std::size_t step) {
    for (int b = 0; b < 3; ++b)
        for (std::size_t j = 0; j < y.blades[b].size(); ++j) {
            diagnostics::LookupLocation location(
                {step * model_.data().dt, std::size_t(b + 1), j + 1, "UA_advance"});
            airfoils_[b][j].advance(y.blades[b][j].alpha, y.blades[b][j].speed, step);
            previous_phi_[b][j] = y.blades[b][j].root_phi;
        }
}
std::array<std::vector<PointLoad>, 3> Rotor::structural_loads(double time, const RotorState &state,
                                                              const RotorOutput &y) const {
    std::array<std::vector<PointLoad>, 3> result;
    for (std::size_t b = 0; b < result.size(); ++b)
        result[b] = structural_loads_for_blade(b, time, state[b], y);
    return result;
}
std::vector<PointLoad> Rotor::structural_loads_for_blade(std::size_t b, double time, const ModalState &state,
                                                         const RotorOutput &y) const {
    std::vector<Motion> motions;
    LoadWorkspace workspace;
    structure_.motions_into(structure_.blade_basis(time, b), state, motions);
    structural_loads_into(b, y, motions, workspace);
    return std::move(workspace.result);
}
void Rotor::structural_loads_into(std::size_t b, const RotorOutput &y, const std::vector<Motion> &motions,
                                  LoadWorkspace &w) const {
    const auto &blade = y.blades.at(b);
    w.source.resize(blade.size());
    w.distributed.resize(blade.size());
    w.destination.resize(structure_.nodes.size() - 2);
    for (std::size_t j = 0; j < blade.size(); ++j) {
        w.source[j] = blade[j].motion.position;
        w.distributed[j] = blade[j].load;
    }
    for (std::size_t j = 1; j + 1 < structure_.nodes.size(); ++j)
        w.destination[j - 1] = motions[j].position;
    load_maps_[b].transfer_into(w.distributed, w.source, w.destination, w.points);
    w.result.assign(structure_.nodes.size(), PointLoad{});
    std::copy(w.points.begin(), w.points.end(), w.result.begin() + 1);
}
Vec3 Rotor::structural_acceleration(std::size_t b, const Matrix3 &basis, const ModalState &state,
                                    const RotorOutput &y, std::vector<Motion> &motions,
                                    LoadWorkspace &workspace) const {
    // Always rebuild for the supplied state, including every Jacobian perturbation.
    structure_.motions_into(basis, state, motions);
    structural_loads_into(b, y, motions, workspace);
    return structure_.acceleration(b, state, motions, workspace.result);
}
} // namespace turbine
