#pragma once
#include "turbine/math.hpp"
namespace turbine {
struct ModalState {
    Vec3 q{}, qd{};
};
struct PointLoad {
    Vec3 force{}, moment{};
};
struct StructuralStation {
    double span = 0, width = 0, mass = 0, twist = 0;
    Vec3 flap{}, edge{}, flap_slope{}, edge_slope{};
    Matrix3 axial{};
};
struct Motion {
    Vec3 position{}, velocity{}, acceleration_bias{}, angular_velocity{};
    Matrix3 orientation{}; // local structural basis in global coordinates
    std::array<Vec3, 3> partial_velocity{}, partial_angular{};
};
class BladeStructure {
  public:
    explicit BladeStructure(const Case &);
    std::vector<StructuralStation> nodes; // root, midpoint quadrature nodes, tip
    Matrix3 stiffness{}, damping{};
    Vec3 modal_mass{};
    double length, hub_radius, omega, initial_azimuth;
    Vec3 hub, shaft;
    std::array<double, 3> cone{}, pitch{}, tip_mass{};
    Matrix3 blade_basis(double time, std::size_t blade) const;
    Motion motion(double time, std::size_t blade, const ModalState &, const StructuralStation &) const;
    Vec3 acceleration(double time, std::size_t blade, const ModalState &,
                      const std::vector<PointLoad> &) const;

  private:
    double gravity_, tilt_, yaw_;
};
} // namespace turbine
