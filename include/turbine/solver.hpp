#pragma once
#include "turbine/rotor.hpp"
namespace turbine {
struct SolverOptions : NewtonOptions {
    // Fixed-base backend modal scales; generic integrators use DofDescriptor scales.
    Vec3 acceleration_scale{1., 1., 1.};
};
class Solver {
  public:
    explicit Solver(const Case &, SolverOptions = {});
    explicit Solver(TurbineModel, SolverOptions = {});
    const SolverOptions &options() const noexcept { return options_; }
    const SolverDiagnostics &diagnostics() const noexcept { return diagnostics_; }
    // Copies share frozen model data; all evolving states and workspaces are independent.
    const Rotor &rotor() const noexcept { return rotor_; }
    const RotorState &state() const noexcept { return state_; }
    const SecondOrderState &generalized_state() const noexcept { return integrator_.state(); }
    const DofLayout &layout() const noexcept { return integrator_.layout(); }
    const RotorOutput &aerodynamic() const noexcept { return aerodynamic_; }
    const std::array<Vec3, 3> &acceleration() const noexcept { return acceleration_; }
    double time() const noexcept { return time_; }
    std::size_t step_number() const noexcept { return step_number_; }
    bool failed() const noexcept { return failed_; }
    void step();
    void reset();
    class Checkpoint {
        friend class Solver;
        std::shared_ptr<const Solver> saved_;
        explicit Checkpoint(const Solver &s) : saved_(std::make_shared<const Solver>(s)) {}

      public:
        Checkpoint(const Checkpoint &) = default;
        Checkpoint &operator=(const Checkpoint &) = default;
    };
    Checkpoint checkpoint() const;
    void restore(const Checkpoint &);

  private:
    Rotor rotor_;
    GeneralizedAlpha integrator_;
    SolverOptions options_;
    SolverDiagnostics diagnostics_;
    RotorState state_{};
    RotorOutput aerodynamic_;
    std::array<Vec3, 3> acceleration_{};
    double time_ = 0;
    std::size_t step_number_ = 0;
    bool failed_ = false;
    RotorWorkspace rotor_workspace_;
    LoadWorkspace load_workspace_;
    double dt_;
    void advance();
};
} // namespace turbine
