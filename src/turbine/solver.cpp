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
    aerodynamic = rotor.evaluate(0, state);
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
    aerodynamic = rotor.evaluate(next_time, predicted);
    std::array<Vec3, 3> current_acc{};
    state = predicted;
    for (int iteration = 0; iteration < 12; ++iteration) {
        const auto loads = rotor.structural_loads(next_time, state, aerodynamic);
        double largest = 0;
        for (int b = 0; b < 3; ++b) {
            const auto f = rotor.structure.acceleration(next_time, b, state[b], loads[b]);
            const Vec3 residual = f - current_acc[b];
            Matrix3 jac{};
            for (int j = 0; j < 3; ++j) {
                const double h = 1e-4;
                auto perturbed = state;
                perturbed[b].q[j] += beta_prime_ * h;
                perturbed[b].qd[j] += gamma_prime_ * h;
                const auto lp = rotor.structural_loads(next_time, perturbed, aerodynamic);
                const auto fp = rotor.structure.acceleration(next_time, b, perturbed[b], lp[b]);
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
