#include "turbine/results.hpp"
#include "acoustic_levels.hpp"
namespace turbine {
OutputLayout::OutputLayout(const Case &c, const AcousticConfiguration &config)
    : parameters(config.parameters),
      prefix(std::filesystem::path(c.acoustic.value("AAOutFile")).filename().string()),
      output_count(c.acoustic.integer("NrOutFile")), blades(config.blades), nodes(config.span.size()),
      first(aeroacoustics::blade_elements(config.span, config.blade_percent).first),
      observers(config.observers.size()), frequencies(parameters.freqlist.size()) {
    for (std::size_t o = 0; o < observers; ++o) {
        const auto observer = "Obs" + std::to_string(o + 1);
        labels[0].push_back(observer);
        for (double f : parameters.freqlist) {
            std::ostringstream value;
            value << f;
            labels[1].push_back(observer + "_Freq" + value.str());
            for (const auto &m : aeroacoustics::mechanism_registry)
                labels[2].push_back(observer + "_Freq" + value.str() + "_" + m.name);
        }
    }
    for (int b = 1; b <= 3; ++b)
        for (std::size_t j = 0; j < nodes; ++j)
            for (std::size_t o = 0; o < observers; ++o)
                labels[3].push_back("Blade" + std::to_string(b) + "_Node" + std::to_string(j + 1) + "_Obs" +
                                    std::to_string(o + 1));
}
AcousticAggregator::AcousticAggregator(OutputLayout layout) : layout_(std::move(layout)) {
    const auto &l = layout_;
    if (l.first >= l.nodes || !l.blades || !l.observers || !l.frequencies)
        throw std::invalid_argument("Invalid aggregation layout");
    result_.power = {std::vector<double>(l.observers), std::vector<double>(l.observers * l.frequencies),
                     std::vector<double>(l.observers * l.frequencies * aeroacoustics::mechanism_count),
                     std::vector<double>(l.blades * l.nodes * l.observers)};
}
const AcousticResult &AcousticAggregator::aggregate(double time, const aeroacoustics::Snapshot &snapshot) {
    const auto no = layout_.observers, nf = layout_.frequencies, first = layout_.first;
    const auto &parameters = layout_.parameters;
    if (snapshot.size() != no)
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
    for (auto &v : values)
        std::fill(v.begin(), v.end(), 0.);
    for (std::size_t o = 0; o < no; ++o)
        for (std::size_t n = 0; n < snapshot[o].size(); ++n) {
            double node_power = 0;
            for (std::size_t m = 0; m < aeroacoustics::mechanism_count; ++m)
                for (std::size_t f = 0; f < nf; ++f) {
                    const double db = snapshot[o][n][m][f];
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
                    mechanisms[(o * nf + f) * aeroacoustics::mechanism_count + m] += power;
                    if (aeroacoustics::mechanism_registry[m].contributes_to_total) {
                        total[o] += power;
                        spectra[o * nf + f] += power;
                        node_power += power;
                    }
                }
            const auto b = n / (layout_.nodes - first), j = n % (layout_.nodes - first) + first;
            nodal[(b * layout_.nodes + j) * no + o] = node_power;
        }

    return result_;
}
} // namespace turbine
