#pragma once
#include "turbine/input.hpp"
namespace turbine {
// A dataset is a measured/computed state, never a roughness-to-noise correction.
struct SurfaceData {
    std::size_t airfoil; // zero based, applies to all stations using this airfoil
    std::string name, provenance, uncertainty_note;
    double alpha_min, alpha_max, re_min, re_max;
    double transition_suction, transition_pressure, roughness_m, erosion_m, relative_uncertainty;
    double te_thickness_m, te_angle_deg;
    std::optional<Airfoil> polar;
    aeroacoustics::BLTable boundary_layer;
    aeroacoustics::PreparedBLTable prepared;
    explicit SurfaceData(const std::filesystem::path &);
    aeroacoustics::BoundaryLayer at(double alpha_deg, double re, double chord) const;
};
class SurfaceSet {
    std::vector<SurfaceData> data_;

  public:
    static std::shared_ptr<const SurfaceSet> read(const std::filesystem::path &);
    const std::vector<SurfaceData> &data() const noexcept { return data_; }
    const SurfaceData *find(std::size_t airfoil) const noexcept;
    Case apply(const Case &) const;
};
} // namespace turbine
