#include "turbine/solver.hpp"
namespace turbine {
Solver::Solver(const Case &c) : Solver(TurbineModel(c)) {}
Solver::Solver(TurbineModel model) : rotor_(std::move(model)), dt_(rotor_.model().data().dt) {
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
void Solver::reset() { *this = Solver(rotor_.model()); }
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
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void Solver::advance() {
    const double next_time = (step_number_ + 1) * dt_;
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
    for (int iteration = 0; iteration < 12; ++iteration) {
        double largest = 0;
        for (int b = 0; b < 3; ++b) {
            const auto f = rotor_.structural_acceleration(b, basis[b], state_[b], aerodynamic_,
                                                          rotor_workspace_.structural, load_workspace_);
            const Vec3 residual = f - current_acc[b];
            Matrix3 jac{};
            for (int j = 0; j < 3; ++j) {
                const double h = 1e-4;
                auto perturbed = state_[b];
                perturbed.q[j] += beta_prime_ * h;
                perturbed.qd[j] += gamma_prime_ * h;
                const auto fp = rotor_.structural_acceleration(b, basis[b], perturbed, aerodynamic_,
                                                               rotor_workspace_.structural, load_workspace_);
                for (int i = 0; i < 3; ++i)
                    jac[i][j] = (i == j ? 1. : 0.) - (fp[i] - f[i]) / h;
            }
            const auto delta = solve3(jac, residual);
            largest = std::max(largest, norm(delta));
            current_acc[b] = current_acc[b] + delta;
            state_[b].q = predicted[b].q + beta_prime_ * current_acc[b];
            state_[b].qd = predicted[b].qd + gamma_prime_ * current_acc[b];
        }
        if (largest < 1e-9)
            break;
        if (iteration == 11)
            throw std::runtime_error("Structural solve failed at time " + std::to_string(next_time));
    }
    for (int b = 0; b < 3; ++b) {
        acceleration_[b] = current_acc[b];
        algorithmic_[b] = a0[b] + ((1 - alpha_f_) / (1 - alpha_m_)) * current_acc[b];
    }
    time_ = next_time;
    ++step_number_;
}
} // namespace turbine
