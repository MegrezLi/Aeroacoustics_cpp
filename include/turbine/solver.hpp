#pragma once
#include "turbine/rotor.hpp"
namespace turbine {
enum class SolverMode { reference, scaled };
struct SolverOptions {
    SolverMode mode = SolverMode::reference;
    int max_iterations = 12;
    double perturbation = 1e-4, correction_tolerance = 1e-9;
    Vec3 acceleration_scale{1., 1., 1.};
    double residual_absolute = 1e-9, residual_relative = 1e-9;
    bool reuse_jacobian = true; // scaled mode only, within one step
};
struct SolverDiagnostics {
    SolverMode mode = SolverMode::reference;
    std::size_t iterations = 0, acceleration_evaluations = 0, jacobian_builds = 0;
    int last_iterations = 0, max_iterations = 0, blade = 0;
    double time = 0, residual = 0, scaled_residual = 0, correction = 0;
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
    SolverOptions options_;
    SolverDiagnostics diagnostics_;
    RotorState state_{};
    RotorOutput aerodynamic_;
    std::array<Vec3, 3> acceleration_{}, algorithmic_{};
    double time_ = 0;
    std::size_t step_number_ = 0;
    bool failed_ = false;
    RotorWorkspace rotor_workspace_;
    LoadWorkspace load_workspace_;
    double dt_, alpha_m_, alpha_f_, beta_, gamma_, beta_prime_, gamma_prime_;
    void advance();
};
} // namespace turbine
