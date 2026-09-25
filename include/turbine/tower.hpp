#pragma once
#include "turbine/input.hpp"
#include <memory>
namespace turbine {
struct TowerStation {
    double z, diameter, cd;
};
// Fixed vertical tower, in the turbine's local coordinates. No tower structural DOFs.
class TowerInfluence {
  public:
    std::string provenance;
    double x = 0, y = 0;
    bool potential = true, shadow = false;
    std::vector<TowerStation> stations;
    void validate() const;
    Vec3 apply(Vec3 wind, const Vec3 &position) const;
    static std::shared_ptr<const TowerInfluence> read(const std::filesystem::path &);
};
} // namespace turbine
