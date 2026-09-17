#pragma once
#include <cstddef>
#include <string>
#include <vector>
namespace turbine {
struct DofDescriptor {
    std::string name, unit;
    double acceleration_scale = 1.;
    std::string output_stem; // Optional legacy CSV stem; otherwise use name.
};
struct DofBlock {
    std::string name;
    std::size_t offset, size;
};
// Each block is one fully coupled system; different blocks must be independent
// during the correction with external/module inputs frozen.
class DofLayout {
    std::vector<DofDescriptor> dofs_;
    std::vector<DofBlock> blocks_;

  public:
    DofLayout(std::vector<DofDescriptor>, std::vector<DofBlock>);
    const auto &dofs() const noexcept { return dofs_; }
    const auto &blocks() const noexcept { return blocks_; }
    std::size_t size() const noexcept { return dofs_.size(); }
};
struct SecondOrderState {
    std::vector<double> q, qd;
};
struct StateView {
    const double *q, *qd;
    std::size_t size;
};
class AccelerationOperator {
  public:
    virtual ~AccelerationOperator() = default;
    // Overwrite all entries. No history advancement here: this may be called
    // repeatedly at trial states. Inputs belong only to the selected block.
    virtual void evaluate(std::size_t block, StateView, double *acceleration) const = 0;
};
enum class SolverMode { reference, scaled };
struct NewtonOptions {
    SolverMode mode = SolverMode::reference;
    int max_iterations = 12;
    double perturbation = 1e-4, correction_tolerance = 1e-9;
    double residual_absolute = 1e-9, residual_relative = 1e-9;
    bool reuse_jacobian = true;
};
struct SolverDiagnostics {
    SolverMode mode = SolverMode::reference;
    std::size_t iterations = 0, acceleration_evaluations = 0, jacobian_builds = 0;
    int last_iterations = 0, max_iterations = 0, blade = 0;
    std::size_t block = 0;
    double time = 0, residual = 0, scaled_residual = 0, correction = 0;
};
// Value-owned state/history/workspace: copies are independent checkpoints.
// A failed correction must be discarded/restored or reset before advancing.
class GeneralizedAlpha {
    struct BlockWorkspace {
        std::vector<double> jac, matrix, rhs, delta, f, fp;
        SecondOrderState trial;
        double previous_residual = 0;
        int age = 0;
        explicit BlockWorkspace(std::size_t);
    };
    DofLayout layout_;
    NewtonOptions options_;
    SolverDiagnostics diagnostics_;
    SecondOrderState state_, predicted_;
    std::vector<double> acceleration_, algorithmic_, a0_, current_;
    std::vector<BlockWorkspace> work_;
    double dt_, alpha_m_, alpha_f_, beta_, gamma_, beta_prime_, gamma_prime_, time_ = 0;
    std::size_t step_ = 0;
    bool pending_ = false, failed_ = false;
    void correct_impl(const AccelerationOperator &);

  public:
    GeneralizedAlpha(DofLayout, double dt, double rho_infinity, NewtonOptions = {});
    const SecondOrderState &predict();
    void correct(const AccelerationOperator &);
    void reset();
    const auto &layout() const noexcept { return layout_; }
    const auto &state() const noexcept { return state_; }
    const auto &acceleration() const noexcept { return acceleration_; }
    const auto &diagnostics() const noexcept { return diagnostics_; }
    double time() const noexcept { return time_; }
    std::size_t step_number() const noexcept { return step_; }
    bool failed() const noexcept { return failed_; }
};
} // namespace turbine
