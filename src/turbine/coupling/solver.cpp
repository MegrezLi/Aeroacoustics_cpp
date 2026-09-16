#include "turbine/solver.hpp"
#include <algorithm>
#include <sstream>
namespace turbine {
Solver::Solver(const Case &c, SolverOptions options) : Solver(TurbineModel(c), options) {}
Solver::Solver(TurbineModel model, SolverOptions options)
    : rotor_(std::move(model)), options_(options), dt_(rotor_.model().data().dt) {
    const auto positive = [](double x) { return std::isfinite(x) && x > 0; };
    if ((options.mode != SolverMode::reference && options.mode != SolverMode::scaled) ||
        options.max_iterations < 1 || !positive(options.perturbation) ||
        !positive(options.correction_tolerance) || !positive(options.residual_absolute) ||
        !std::isfinite(options.residual_relative) || options.residual_relative < 0 ||
        !std::all_of(options.acceleration_scale.begin(), options.acceleration_scale.end(), positive))
        throw std::invalid_argument("Invalid structural solver options");
    diagnostics_.mode = options_.mode;
    const auto &c = rotor_.model().data();
    if (c.primary.integer("ModCoupling") != 3 || c.primary.integer("NumCrctn") != 0)
        throw std::runtime_error("Expected ModCoupling=3 and NumCrctn=0");
    const double rho = c.primary.number("RhoInf");
    alpha_m_ = (2 * rho - 1) / (rho + 1);
    alpha_f_ = rho / (rho + 1);
    gamma_ = .5 - alpha_m_ + alpha_f_;
    beta_ = std::pow(1 - alpha_m_ + alpha_f_, 2) / 4;
    beta_prime_ = dt_ * dt_ * beta_ * (1 - alpha_f_) / (1 - alpha_m_);
    gamma_prime_ = dt_ * gamma_ * (1 - alpha_f_) / (1 - alpha_m_);
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
    diagnostics_.time = next_time;
    diagnostics_.last_iterations = diagnostics_.blade = 0;
    diagnostics_.residual = diagnostics_.scaled_residual = diagnostics_.correction = 0;
    RotorState predicted = state_;
    std::array<Vec3, 3> a0{};
    for (int b = 0; b < 3; ++b) {
        a0[b] = (alpha_f_ * acceleration_[b] - alpha_m_ * algorithmic_[b]) / (1 - alpha_m_);
        predicted[b].qd = state_[b].qd + dt_ * (1 - gamma_) * algorithmic_[b] + gamma_ * dt_ * a0[b];
        predicted[b].q = state_[b].q + dt_ * state_[b].qd + dt_ * dt_ * (.5 - beta_) * algorithmic_[b] +
                         beta_ * dt_ * dt_ * a0[b];
    }
    // AeroDyn is an Option 2 module in FAST_Solver: it advances once from
    // the predicted structural motion before the structural Newton solve.
    rotor_.advance_airfoils(aerodynamic_, step_number_);
    rotor_.evaluate_into(next_time, predicted, aerodynamic_, rotor_workspace_);
    std::array<Matrix3, 3> basis;
    for (int b = 0; b < 3; ++b)
        basis[b] = rotor_.structure().blade_basis(next_time, b);
    std::array<Vec3, 3> current_acc{};
    state_ = predicted;
    std::array<Matrix3, 3> jacobians{};
    std::array<double, 3> previous_residual{};
    std::array<int, 3> age{};
    const bool scaled = options_.mode == SolverMode::scaled;
    auto finite = [](const Vec3 &v) {
        for (double x : v)
            if (!std::isfinite(x))
                throw std::runtime_error("Non-finite structural value");
    };
    auto evaluate = [&](int b, const ModalState &state) {
        ++diagnostics_.acceleration_evaluations;
        finite(state.q);
        finite(state.qd);
        auto f = rotor_.structural_acceleration(b, basis[b], state, aerodynamic_, rotor_workspace_.structural,
                                                load_workspace_);
        finite(f);
        return f;
    };
    for (int iteration = 0; iteration < options_.max_iterations; ++iteration) {
        ++diagnostics_.iterations;
        diagnostics_.last_iterations = iteration + 1;
        diagnostics_.max_iterations = std::max(diagnostics_.max_iterations, iteration + 1);
        double largest = 0, worst_residual = 0, worst_scaled = 0;
        for (int b = 0; b < 3; ++b) {
            diagnostics_.blade = b + 1;
            const auto f = evaluate(b, state_[b]);
            const Vec3 residual = f - current_acc[b];
            finite(residual);
            double scaled_residual = 0;
            for (int i = 0; i < 3; ++i) {
                const double scale =
                    std::max({options_.acceleration_scale[i], std::abs(f[i]), std::abs(current_acc[b][i])});
                const double tolerance = options_.residual_absolute + options_.residual_relative * scale;
                if (!std::isfinite(tolerance) || tolerance <= 0)
                    throw std::runtime_error("Invalid scaled residual tolerance");
                scaled_residual = std::max(scaled_residual, std::abs(residual[i]) / tolerance);
            }
            diagnostics_.residual = norm(residual);
            diagnostics_.scaled_residual = scaled_residual;
            worst_residual = std::max(worst_residual, diagnostics_.residual);
            worst_scaled = std::max(worst_scaled, scaled_residual);
            auto &jac = jacobians[b];
            // Reuse only within this step. Stagnation/worsening or two reuses forces a rebuild.
            if (!scaled || !options_.reuse_jacobian || iteration == 0 || age[b] >= 2 ||
                scaled_residual >= .9 * previous_residual[b]) {
                ++diagnostics_.jacobian_builds;
                age[b] = 0;
                for (int j = 0; j < 3; ++j) {
                    const double h =
                        options_.perturbation *
                        (scaled ? std::max(options_.acceleration_scale[j], std::abs(current_acc[b][j])) : 1.);
                    if (!std::isfinite(h) || h <= 0)
                        throw std::runtime_error("Invalid Jacobian perturbation");
                    auto perturbed = state_[b];
                    perturbed.q[j] += beta_prime_ * h;
                    perturbed.qd[j] += gamma_prime_ * h;
                    const auto fp = evaluate(b, perturbed);
                    for (int i = 0; i < 3; ++i)
                        jac[i][j] = (i == j ? 1. : 0.) - (fp[i] - f[i]) / h;
                }
                for (const auto &row : jac)
                    finite(row);
            } else
                ++age[b];
            previous_residual[b] = scaled_residual;
            const auto delta = solve3(jac, residual);
            finite(delta);
            const double correction = norm(delta);
            if (!std::isfinite(correction) || !std::isfinite(diagnostics_.residual) ||
                !std::isfinite(scaled_residual))
                throw std::runtime_error("Non-finite structural norm");
            diagnostics_.correction = correction;
            largest = std::max(largest, correction);
            current_acc[b] = current_acc[b] + delta;
            state_[b].q = predicted[b].q + beta_prime_ * current_acc[b];
            state_[b].qd = predicted[b].qd + gamma_prime_ * current_acc[b];
            finite(state_[b].q);
            finite(state_[b].qd);
        }
        diagnostics_.residual = worst_residual;
        diagnostics_.scaled_residual = worst_scaled;
        diagnostics_.correction = largest;
        if (largest < options_.correction_tolerance && (!scaled || worst_scaled <= 1.))
            break;
        if (iteration + 1 == options_.max_iterations)
            throw std::runtime_error("Structural iteration limit exceeded (reported norms are blade maxima)");
    }
    for (int b = 0; b < 3; ++b) {
        acceleration_[b] = current_acc[b];
        algorithmic_[b] = a0[b] + ((1 - alpha_f_) / (1 - alpha_m_)) * current_acc[b];
    }
    time_ = next_time;
    ++step_number_;
}
} // namespace turbine
