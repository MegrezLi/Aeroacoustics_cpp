#include "source_models.hpp"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace aeroacoustics {
struct AcousticWorkspace::Impl {
    Parameters parameters;
    SourceSelection models;
    Spectrum weighting;
    detail::TnoWorkspace integration;
    std::vector<detail::PreparedSection> sources;
    Snapshot output;

    explicit Impl(Parameters p) : parameters(std::move(p)), models(parameters) {
        if (parameters.tbltemod == 2)
            parameters.x_blmethod = 2;
        if (parameters.aweighting)
            weighting = a_weighting(parameters.freqlist);
    }
};
AcousticWorkspace::AcousticWorkspace(Parameters p) : impl_(std::make_unique<Impl>(std::move(p))) {}
AcousticWorkspace::~AcousticWorkspace() = default;
AcousticWorkspace::AcousticWorkspace(const AcousticWorkspace &other)
    : impl_(std::make_unique<Impl>(*other.impl_)) {}
AcousticWorkspace &AcousticWorkspace::operator=(const AcousticWorkspace &other) {
    if (this != &other)
        impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}
AcousticWorkspace::AcousticWorkspace(AcousticWorkspace &&) noexcept = default;
AcousticWorkspace &AcousticWorkspace::operator=(AcousticWorkspace &&) noexcept = default;

const Snapshot &AcousticWorkspace::evaluate(const std::vector<Node> &nodes,
                                            const std::vector<Vec3> &observers) {
    evaluate_blocks(
        nodes, observers, [](std::size_t, const Snapshot &) {}, std::max<std::size_t>(1, observers.size()));
    return impl_->output;
}
void AcousticWorkspace::evaluate_blocks(const std::vector<Node> &nodes, const std::vector<Vec3> &observers,
                                        const ObserverBlockCallback &consume, std::size_t block_size) {
    if (!consume || !block_size)
        throw std::invalid_argument("Invalid observer block request");
    auto &work = *impl_;
    if (observers.empty()) {
        work.output.clear();
        return;
    }
    work.sources.resize(nodes.size());
    // Boundary layers and all source spectral shapes depend on the section,
    // not on the observer. Prepare them exactly once for this snapshot.
    for (std::size_t n = 0; n < nodes.size(); ++n)
        detail::prepare_section(work.parameters, work.models, nodes[n].section, work.integration,
                                work.sources[n]);
    // Shrink if this workspace previously served the full-snapshot API.
    if (work.output.size() > block_size)
        work.output.resize(block_size);
    for (std::size_t first = 0; first < observers.size();) {
        const auto count = std::min(block_size, observers.size() - first);
        work.output.resize(count);
        for (std::size_t local = 0; local < count; ++local) {
            const auto o = first + local;
            auto &output = work.output[local];
            output.resize(nodes.size());
            for (std::size_t n = 0; n < nodes.size(); ++n) {
                const auto &node = nodes[n];
                const auto geometry = observe(observers[o], node.aero_center, node.global_to_local,
                                              node.section.chord, node.airfoil_reference);
                try {
                    detail::emit_section(work.parameters, work.models, work.sources[n], geometry.first,
                                         geometry.second, work.weighting, output[n]);
                } catch (const std::exception &e) {
                    throw std::runtime_error("observer=" + std::to_string(o + 1) +
                                             " blade=" + std::to_string(node.blade_number) +
                                             " node=" + std::to_string(node.node_number) +
                                             " selected_node=" + std::to_string(n + 1) + ": " + e.what());
                }
            }
        }
        consume(first, work.output);
        first += count;
    }
}
} // namespace aeroacoustics
