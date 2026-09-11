#include "source_models.hpp"
#include <utility>
#include <stdexcept>

namespace aeroacoustics {
struct AcousticWorkspace::Impl {
    Parameters parameters;
    Spectrum weighting;
    detail::TnoWorkspace integration;
    std::vector<detail::PreparedSection> sources;
    Snapshot output;

    explicit Impl(Parameters p) : parameters(std::move(p)) {
        validate(parameters);
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
    auto &work = *impl_;
    work.sources.resize(nodes.size());
    // Boundary layers and all source spectral shapes depend on the section,
    // not on the observer. Prepare them exactly once for this snapshot.
    for (std::size_t n = 0; n < nodes.size(); ++n)
        detail::prepare_section(work.parameters, nodes[n].section, work.integration, work.sources[n]);
    work.output.resize(observers.size());
    for (std::size_t o = 0; o < observers.size(); ++o) {
        auto &output = work.output[o];
        output.resize(nodes.size());
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            const auto &node = nodes[n];
            const auto geometry = observe(observers[o], node.aero_center, node.global_to_local,
                                          node.section.chord, node.airfoil_reference);
            try {
                detail::emit_section(work.parameters, work.sources[n], geometry.first, geometry.second,
                                     work.weighting, output[n]);
            } catch (const std::exception &e) {
                throw std::runtime_error("observer=" + std::to_string(o + 1) + " blade=" +
                    std::to_string(node.blade_number) + " node=" + std::to_string(node.node_number) +
                    " selected_node=" + std::to_string(n + 1) + ": " + e.what());
            }
        }
    }
    return work.output;
}
} // namespace aeroacoustics
