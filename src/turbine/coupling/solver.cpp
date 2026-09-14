#include "turbine/solver.hpp"
namespace turbine {
Solver::Solver(const Case &c) : rotor(c), dt_(c.dt) {
    if (c.primary.integer("ModCoupling") != 3 || c.primary.integer("NumCrctn") != 0)
        throw std::runtime_error("Expected ModCoupling=3 and NumCrctn=0");
    const double rho = c.primary.number("RhoInf");
    alpha_m_ = (2 * rho - 1) / (rho + 1);
    alpha_f_ = rho / (rho + 1);
    gamma_ = .5 - alpha_m_ + alpha_f_;
    beta_ = std::pow(1 - alpha_m_ + alpha_f_, 2) / 4;
    beta_prime_ = dt_ * dt_ * beta_ * (1 - alpha_f_) / (1 - alpha_m_);
    gamma_prime_ = dt_ * gamma_ * (1 - alpha_f_) / (1 - alpha_m_);
    rotor.evaluate_into(0, state, aerodynamic, rotor_workspace_);
    // FAST_Solver initializes its acceleration arrays to zero; Step0 solves
    // module inputs without populating these arrays for ElastoDyn.
}
void Solver::step() {
    const double next_time = (step_number + 1) * dt_;
    RotorState predicted = state;
    std::array<Vec3, 3> a0{};
    for (int b = 0; b < 3; ++b) {
        a0[b] = (alpha_f_ * acceleration[b] - alpha_m_ * algorithmic[b]) / (1 - alpha_m_);
        predicted[b].qd = state[b].qd + dt_ * (1 - gamma_) * algorithmic[b] + gamma_ * dt_ * a0[b];
        predicted[b].q = state[b].q + dt_ * state[b].qd + dt_ * dt_ * (.5 - beta_) * algorithmic[b] +
                         beta_ * dt_ * dt_ * a0[b];
    }
    // AeroDyn is an Option 2 module in FAST_Solver: it advances once from
    // the predicted structural motion before the structural Newton solve.
    rotor.advance_airfoils(aerodynamic, step_number);
    rotor.evaluate_into(next_time, predicted, aerodynamic, rotor_workspace_);
    std::array<Matrix3, 3> basis;
    for (int b = 0; b < 3; ++b)
        basis[b] = rotor.structure.blade_basis(next_time, b);
    std::array<Vec3, 3> current_acc{};
    state = predicted;
    for (int iteration = 0; iteration < 12; ++iteration) {
        double largest = 0;
        for (int b = 0; b < 3; ++b) {
            const auto f = rotor.structural_acceleration(b, basis[b], state[b], aerodynamic,
                                                         rotor_workspace_.structural, load_workspace_);
            const Vec3 residual = f - current_acc[b];
            Matrix3 jac{};
            for (int j = 0; j < 3; ++j) {
                const double h = 1e-4;
                auto perturbed = state[b];
                perturbed.q[j] += beta_prime_ * h;
                perturbed.qd[j] += gamma_prime_ * h;
                const auto fp = rotor.structural_acceleration(b, basis[b], perturbed, aerodynamic,
                                                              rotor_workspace_.structural, load_workspace_);
                for (int i = 0; i < 3; ++i)
                    jac[i][j] = (i == j ? 1. : 0.) - (fp[i] - f[i]) / h;
            }
            const auto delta = solve3(jac, residual);
            largest = std::max(largest, norm(delta));
            current_acc[b] = current_acc[b] + delta;
            state[b].q = predicted[b].q + beta_prime_ * current_acc[b];
            state[b].qd = predicted[b].qd + gamma_prime_ * current_acc[b];
        }
        if (largest < 1e-9)
            break;
        if (iteration == 11)
            throw std::runtime_error("Structural solve failed at time " + std::to_string(next_time));
    }
    for (int b = 0; b < 3; ++b) {
        acceleration[b] = current_acc[b];
        algorithmic[b] = a0[b] + ((1 - alpha_f_) / (1 - alpha_m_)) * current_acc[b];
    }
    time = next_time;
    ++step_number;
}
} // namespace turbine
