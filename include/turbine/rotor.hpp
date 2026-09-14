#pragma once
#include "turbine/bem.hpp"
#include "turbine/mesh.hpp"
#include "turbine/unsteady.hpp"
namespace turbine {
using RotorState = std::array<ModalState, 3>;
struct AeroStation {
    Motion motion;
    Vec3 wind;
    Matrix3 annulus;
    BEMInput bem;
    double root_phi = 0, phi = 0, alpha = 0, speed = 0, axial = 0, tangential = 0;
    Coefficients coefficients;
    PointLoad load;
};
struct RotorOutput {
    std::array<std::vector<AeroStation>, 3> blades;
    Vec3 average_velocity{};
    double skew = 0;
};
// Owned by the caller/Solver, never shared as hidden mutable Rotor state.
struct RotorWorkspace {
    std::vector<Motion> structural, aerodynamic;
};
struct LoadWorkspace {
    std::vector<Vec3> source, destination;
    std::vector<PointLoad> distributed, points, result;
};
class Rotor {
  public:
    explicit Rotor(const Case &);
    BladeStructure structure;
    RotorOutput evaluate(double time, const RotorState &) const;
    void evaluate_into(double time, const RotorState &, RotorOutput &, RotorWorkspace &) const;
    // Reuses motions for both load mapping and acceleration in one state evaluation.
    Vec3 structural_acceleration(std::size_t blade, const Matrix3 &basis, const ModalState &,
                                 const RotorOutput &, std::vector<Motion> &, LoadWorkspace &) const;
    void advance_airfoils(const RotorOutput &, std::size_t step);
    std::array<std::vector<PointLoad>, 3> structural_loads(double time, const RotorState &,
                                                           const RotorOutput &) const;
    std::vector<PointLoad> structural_loads_for_blade(std::size_t blade, double time, const ModalState &,
                                                      const RotorOutput &) const;

  private:
    void structural_loads_into(std::size_t blade, const RotorOutput &, const std::vector<Motion> &,
                               LoadWorkspace &) const;
    const Case *case_;
    BEMOptions options_;
    struct SkewOptions {
        bool redistribute;
        double factor;
    } skew_;
    std::vector<double> tip_constant_, hub_constant_;
    std::vector<MotionMap> motion_maps_;
    std::vector<LoadMap> load_maps_;
    std::array<std::vector<UnsteadyAirfoil>, 3> airfoils_;
    std::array<std::vector<double>, 3> previous_phi_;
};
} // namespace turbine
