#include "turbine/results.hpp"
#include "acoustic_levels.hpp"
namespace turbine {
OutputLayout::OutputLayout(const Case &c, const AcousticConfiguration &config)
    : parameters(config.parameters),
      prefix(std::filesystem::path(c.acoustic.value("AAOutFile")).filename().string()),
      output_count(c.acoustic.integer("NrOutFile")), blades(config.blades), nodes(config.span.size()),
      first(aeroacoustics::blade_elements(config.span, config.blade_percent).first),
      observers(config.observers.size()), frequencies(parameters.freqlist.size()) {
    const auto configuration = configure_modules(c);
    module_profile = configuration.profile;
    const auto layout = structural_dof_layout(configuration);
    dofs = layout.dofs();
    coupling_blocks = layout.blocks();
    for (const auto &block : layout.blocks())
        for (bool velocity : {false, true})
            for (std::size_t i = block.offset; i < block.offset + block.size; ++i)
                dynamics.push_back({i, velocity,
                                    std::string(velocity ? "qd" : "q") +
                                        (dofs[i].output_stem.empty() ? dofs[i].name : dofs[i].output_stem)});
    acoustic_metadata =
        std::make_shared<const aeroacoustics::AcousticMetadata>(aeroacoustics::AcousticMetadata{
            aeroacoustics::FrequencyBands::openfast_reference(parameters.freqlist),
            parameters.aweighting ? aeroacoustics::Weighting::a : aeroacoustics::Weighting::unweighted});
    if (output_count < 1 || output_count > 4)
        throw std::invalid_argument("NrOutFile must be between 1 and 4");
    for (std::size_t o = 0; o < observers; ++o) {
        const auto observer = "Obs" + std::to_string(o + 1);
        labels[0].push_back(observer);
        if (output_count >= 2)
            for (double f : parameters.freqlist) {
                std::ostringstream value;
                value << f;
                labels[1].push_back(observer + "_Freq" + value.str());
                if (output_count >= 3)
                    for (const auto &m : aeroacoustics::mechanism_registry)
                        labels[2].push_back(observer + "_Freq" + value.str() + "_" + m.name);
            }
    }
    if (output_count >= 4)
        for (std::size_t b = 1; b <= blades; ++b)
            for (std::size_t j = 0; j < nodes; ++j)
                for (std::size_t o = 0; o < observers; ++o)
                    labels[3].push_back("Blade" + std::to_string(b) + "_Node" + std::to_string(j + 1) +
                                        "_Obs" + std::to_string(o + 1));
}
AcousticAggregator::AcousticAggregator(OutputLayout layout) : layout_(std::move(layout)) {
    const auto &l = layout_;
    if (l.frequencies != l.parameters.freqlist.size())
        throw std::invalid_argument("Aggregation frequency metadata mismatch");
    result_.metadata =
        std::make_shared<const aeroacoustics::AcousticMetadata>(aeroacoustics::AcousticMetadata{
            aeroacoustics::FrequencyBands::openfast_reference(l.parameters.freqlist),
            l.parameters.aweighting ? aeroacoustics::Weighting::a : aeroacoustics::Weighting::unweighted});
    if (l.first >= l.nodes || !l.blades || !l.observers || !l.frequencies)
        throw std::invalid_argument("Invalid aggregation layout");
    if (l.output_count < 1 || l.output_count > 4)
        throw std::invalid_argument("Invalid aggregation output count");
    const std::array<std::size_t, 4> sizes{l.observers, l.observers * l.frequencies,
                                           l.observers * l.frequencies * aeroacoustics::mechanism_count,
                                           l.blades * l.nodes * l.observers};
    for (int k = 0; k < l.output_count; ++k)
        result_.power[k].resize(sizes[k]);
}
const AcousticResult &AcousticAggregator::aggregate(double time, const aeroacoustics::Snapshot &snapshot) {
    begin();
    append(time, 0, snapshot);
    return finish();
}
void AcousticAggregator::begin() {
    next_observer_ = 0;
    for (auto &v : result_.power)
        std::fill(v.begin(), v.end(), 0.);
}
const AcousticResult &AcousticAggregator::finish() const {
    if (next_observer_ != layout_.observers)
        throw std::logic_error("Incomplete observer aggregation");
    return result_;
}
void AcousticAggregator::append(double time, std::size_t offset, const aeroacoustics::Snapshot &snapshot) {
    const auto no = layout_.observers, nf = layout_.frequencies, first = layout_.first;
    const auto &parameters = layout_.parameters;
    if (offset != next_observer_ || offset > no || snapshot.empty() || snapshot.size() > no - offset)
        throw std::invalid_argument("Aggregation observer count mismatch");
    for (const auto &observer : snapshot) {
        if (observer.size() != layout_.blades * (layout_.nodes - first))
            throw std::invalid_argument("Aggregation node count mismatch");
        for (const auto &node : observer)
            for (const auto &mechanism : node)
                if (mechanism.size() != nf)
                    throw std::invalid_argument("Aggregation frequency count mismatch");
    }
    auto &values = result_.power;
    auto &total = values[0], &spectra = values[1], &mechanisms = values[2], &nodal = values[3];
    for (std::size_t local = 0; local < snapshot.size(); ++local) {
        const auto o = offset + local;
        for (std::size_t n = 0; n < snapshot[local].size(); ++n) {
            double node_power = 0;
            for (std::size_t m = 0; m < aeroacoustics::mechanism_count; ++m)
                for (std::size_t f = 0; f < nf; ++f) {
                    const double db = snapshot[local][n][m][f];
                    double power;
                    try {
                        power = aeroacoustics::relative_power(db);
                    } catch (const std::exception &e) {
                        throw std::runtime_error(std::string(e.what()) + " time=" + std::to_string(time) +
                                                 " observer=" + std::to_string(o + 1) + " blade=" +
                                                 std::to_string(n / (layout_.nodes - first) + 1) + " node=" +
                                                 std::to_string(n % (layout_.nodes - first) + first + 1) +
                                                 aeroacoustics::spectrum_context(parameters, m, f));
                    }
                    if (layout_.output_count >= 3)
                        mechanisms[(o * nf + f) * aeroacoustics::mechanism_count + m] += power;
                    if (aeroacoustics::mechanism_registry[m].contributes_to_total) {
                        total[o] += power;
                        if (layout_.output_count >= 2)
                            spectra[o * nf + f] += power;
                        if (layout_.output_count >= 4)
                            node_power += power;
                    }
                }
            const auto b = n / (layout_.nodes - first), j = n % (layout_.nodes - first) + first;
            if (layout_.output_count >= 4)
                nodal[(b * layout_.nodes + j) * no + o] = node_power;
        }
    }
    next_observer_ += snapshot.size();
}
} // namespace turbine
