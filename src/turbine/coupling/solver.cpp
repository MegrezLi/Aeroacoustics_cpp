#include "turbine/solver.hpp"
#include "turbine/structural_adapter.hpp"
#include <algorithm>
#include <sstream>
namespace turbine {
Solver::Solver(const Case &c, SolverOptions options) : Solver(TurbineModel(c), options) {}
Solver::Solver(TurbineModel model, SolverOptions options)
    : rotor_(std::move(model)),
      integrator_(structural_dof_layout(rotor_.model().configuration(), options.acceleration_scale),
                  rotor_.model().data().dt, rotor_.model().data().primary.number("RhoInf"), options),
      options_(options), dt_(rotor_.model().data().dt) {
    diagnostics_ = integrator_.diagnostics();
    rotor_.evaluate_into(0, state_, aerodynamic_, rotor_workspace_);
    // FAST_Solver initializes its acceleration arrays to zero; Step0 solves
    // module inputs without populating these arrays for ElastoDyn.
}
void Solver::reset() { *this = Solver(rotor_.model(), options_); }
Solver::Checkpoint Solver::checkpoint() const {
    if (failed_)
        throw std::logic_error("Cannot checkpoint a failed solver");
    return Checkpoint(*this);
}
void Solver::restore(const Checkpoint &snapshot) {
    if (!rotor_.model().same_model(snapshot.saved_->rotor_.model()))
        throw std::invalid_argument("Checkpoint belongs to a different turbine model");
    Solver restored = *snapshot.saved_;
    *this = std::move(restored);
}
void Solver::step() {
    if (failed_)
        throw std::logic_error("Solver failed; reset or restore before stepping");
    try {
        advance();
    } catch (const std::exception &e) {
        failed_ = true;
        diagnostics_ = integrator_.diagnostics();
        diagnostics_.blade = static_cast<int>(diagnostics_.block);
        std::ostringstream message;
        message.precision(17);
        message << "Structural/coupling solve time=" << (step_number_ + 1) * dt_
                << " blade=" << diagnostics_.blade << " iteration=" << diagnostics_.last_iterations
                << " residual=" << diagnostics_.residual
                << " scaled_residual=" << diagnostics_.scaled_residual
                << " correction=" << diagnostics_.correction << ": " << e.what();
        throw std::runtime_error(message.str());
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void Solver::advance() {
    const double next_time = (step_number_ + 1) * dt_;
    const auto predicted = rotor_state(integrator_.predict());
    // These module histories advance exactly once. Newton trial evaluations
    // only read the frozen aerodynamic output and do not advance UA or BEM.
    rotor_.advance_airfoils(aerodynamic_, step_number_);
    rotor_.evaluate_into(next_time, predicted, aerodynamic_, rotor_workspace_);
    FixedBaseAcceleration structural(rotor_, aerodynamic_, next_time, rotor_workspace_, load_workspace_);
    integrator_.correct(structural);
    state_ = rotor_state(integrator_.state());
    acceleration_ = rotor_acceleration(integrator_.acceleration());
    diagnostics_ = integrator_.diagnostics();
    diagnostics_.blade = static_cast<int>(diagnostics_.block);
    time_ = integrator_.time();
    step_number_ = integrator_.step_number();
}
} // namespace turbine
