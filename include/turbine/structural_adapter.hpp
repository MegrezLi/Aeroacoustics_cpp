#pragma once
#include "turbine/rotor.hpp"
namespace turbine {
// Short-lived evaluation adapter. All aerodynamic inputs are frozen for the
// current structural correction; scratch buffers are owned by the Solver.
class FixedBaseAcceleration final : public AccelerationOperator {
    const Rotor &rotor_;
    const RotorOutput &aerodynamic_;
    RotorWorkspace &motions_;
    LoadWorkspace &loads_;
    std::array<Matrix3, FixedBaseBladeBackend::blades> basis_;

  public:
    FixedBaseAcceleration(const Rotor &, const RotorOutput &, double time, RotorWorkspace &, LoadWorkspace &);
    void evaluate(std::size_t block, StateView, double *) const override;
};
RotorState rotor_state(const SecondOrderState &);
std::array<Vec3, FixedBaseBladeBackend::blades> rotor_acceleration(const std::vector<double> &);
} // namespace turbine
