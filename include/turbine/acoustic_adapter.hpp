#pragma once
#include "turbine/rotor.hpp"
#include "turbine/surface.hpp"
namespace turbine {
struct AcousticConfiguration {
    aeroacoustics::Parameters parameters;
    aeroacoustics::Spectrum span;
    std::vector<aeroacoustics::Vec3> observers;
    double sample_interval, start, blade_percent, hub_height;
    int ti_method;
    std::size_t blades = 3;
    explicit AcousticConfiguration(const Case &);
    aeroacoustics::AcousticDriver make_driver() const;
};
class AcousticInputAdapter {
    std::vector<std::vector<aeroacoustics::Node>> nodes_;
    std::vector<std::optional<aeroacoustics::PreparedBLTable>> tables_;
    double viscosity_;
    std::shared_ptr<const SurfaceSet> surfaces_;
    std::vector<std::size_t> airfoil_ids_;

  public:
    AcousticInputAdapter(const Case &, const AcousticConfiguration &, std::shared_ptr<const SurfaceSet> = {});
    // Borrowed until the next update. No model pointers are retained.
    const std::vector<std::vector<aeroacoustics::Node>> &update(const RotorOutput &, double time,
                                                                bool sampling, std::size_t first_node);
};
} // namespace turbine
