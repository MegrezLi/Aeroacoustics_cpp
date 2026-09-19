#pragma once
#include "turbine/math.hpp"
#include <memory>
namespace turbine {
class WindField {
  public:
    virtual ~WindField() = default;
    virtual Vec3 at(double time, const Vec3 &position) const = 0;
    virtual const char *description() const noexcept = 0;
};
// Global Cartesian components [m/s], axes time [s], x/y/z [m].
// Values are time-major, then x, y, z (z fastest). No extrapolation.
class GridWind final : public WindField {
    std::array<std::vector<double>, 4> axes_;
    std::vector<Vec3> values_;

  public:
    GridWind(std::array<std::vector<double>, 4>, std::vector<Vec3>);
    static std::shared_ptr<const GridWind> read(const std::filesystem::path &);
    Vec3 at(double, const Vec3 &) const override;
    const char *description() const noexcept override {
        return "global Cartesian time/x/y/z grid; multilinear interpolation";
    }
};
} // namespace turbine
