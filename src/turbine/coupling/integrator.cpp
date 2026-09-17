#include "turbine/integrator.hpp"
#include "turbine/math.hpp"
#include <limits>
#include <set>
namespace turbine {
namespace {
bool positive(double x) { return std::isfinite(x) && x > 0; }
void finite(double x) {
    if (!std::isfinite(x))
        throw std::runtime_error("Non-finite structural value");
}
double vector_norm(const std::vector<double> &v) {
    double value = v.front() * v.front();
    for (std::size_t i = 1; i < v.size(); ++i)
        value += v[i] * v[i];
    return std::sqrt(value);
}
void solve(std::vector<double> &matrix, std::vector<double> &rhs, std::vector<double> &delta) {
    const auto n = rhs.size();
    if (n == 3) {
        Matrix3 a{};
        Vec3 b{};
        for (std::size_t i = 0; i < 3; ++i) {
            b[i] = rhs[i];
            for (std::size_t j = 0; j < 3; ++j)
                a[i][j] = matrix[i * 3 + j];
        }
        const auto x = solve3(a, b);
        std::copy(x.begin(), x.end(), delta.begin());
        return;
    }
    double scale = 0;
    for (double x : matrix) {
        finite(x);
        scale = std::max(scale, std::abs(x));
    }
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t pivot = i;
        for (std::size_t j = i + 1; j < n; ++j)
            if (std::abs(matrix[j * n + i]) > std::abs(matrix[pivot * n + i]))
                pivot = j;
        if (std::abs(matrix[pivot * n + i]) <= std::numeric_limits<double>::epsilon() * scale * n)
            throw std::runtime_error("Singular structural Jacobian");
        for (std::size_t k = 0; k < n; ++k)
            std::swap(matrix[i * n + k], matrix[pivot * n + k]);
        std::swap(rhs[i], rhs[pivot]);
        for (std::size_t j = i + 1; j < n; ++j) {
            const double f = matrix[j * n + i] / matrix[i * n + i];
            for (std::size_t k = i; k < n; ++k)
                matrix[j * n + k] -= f * matrix[i * n + k];
            rhs[j] -= f * rhs[i];
        }
    }
    for (std::size_t ii = n; ii > 0; --ii) {
        const auto i = ii - 1;
        double value = rhs[i];
        for (std::size_t j = i + 1; j < n; ++j)
            value -= matrix[i * n + j] * delta[j];
        delta[i] = value / matrix[i * n + i];
    }
}
} // namespace
DofLayout::DofLayout(std::vector<DofDescriptor> dofs, std::vector<DofBlock> blocks)
    : dofs_(std::move(dofs)), blocks_(std::move(blocks)) {
    if (dofs_.empty() || blocks_.empty())
        throw std::invalid_argument("Empty DOF layout");
    std::set<std::string> names;
    for (const auto &d : dofs_)
        if (d.name.empty() || d.unit.empty() || !positive(d.acceleration_scale) ||
            !names.insert(d.name).second)
            throw std::invalid_argument("Invalid or duplicate DOF descriptor");
    std::size_t offset = 0;
    names.clear();
    for (const auto &b : blocks_) {
        if (b.name.empty() || !names.insert(b.name).second || b.offset != offset || !b.size ||
            b.size > dofs_.size() - offset)
            throw std::invalid_argument("DOF blocks must partition the full layout in order");
        offset += b.size;
    }
    if (offset != dofs_.size())
        throw std::invalid_argument("Incomplete DOF block coverage");
}
GeneralizedAlpha::BlockWorkspace::BlockWorkspace(std::size_t n)
    : jac(n * n), matrix(n * n), rhs(n), delta(n), f(n), fp(n),
      trial{std::vector<double>(n), std::vector<double>(n)} {}
GeneralizedAlpha::GeneralizedAlpha(DofLayout layout, double dt, double rho, NewtonOptions options)
    : layout_(std::move(layout)), options_(options),
      state_{std::vector<double>(layout_.size()), std::vector<double>(layout_.size())}, predicted_(state_),
      acceleration_(layout_.size()), algorithmic_(layout_.size()), a0_(layout_.size()),
      current_(layout_.size()), dt_(dt) {
    if (!positive(dt) || !std::isfinite(rho) || rho < 0 || rho > 1 ||
        (options.mode != SolverMode::reference && options.mode != SolverMode::scaled) ||
        options.max_iterations < 1 || !positive(options.perturbation) ||
        !positive(options.correction_tolerance) || !positive(options.residual_absolute) ||
        !std::isfinite(options.residual_relative) || options.residual_relative < 0)
        throw std::invalid_argument("Invalid generalized-alpha/Newton options");
    diagnostics_.mode = options.mode;
    alpha_m_ = (2 * rho - 1) / (rho + 1);
    alpha_f_ = rho / (rho + 1);
    gamma_ = .5 - alpha_m_ + alpha_f_;
    beta_ = std::pow(1 - alpha_m_ + alpha_f_, 2) / 4;
    beta_prime_ = dt_ * dt_ * beta_ * (1 - alpha_f_) / (1 - alpha_m_);
    gamma_prime_ = dt_ * gamma_ * (1 - alpha_f_) / (1 - alpha_m_);
    if (!positive(beta_prime_) || !positive(gamma_prime_))
        throw std::invalid_argument("Time step produces unrepresentable integration coefficients");
    for (const auto &b : layout_.blocks())
        work_.emplace_back(b.size);
}
const SecondOrderState &GeneralizedAlpha::predict() {
    if (failed_ || pending_)
        throw std::logic_error("Integrator requires reset/restore or completion of pending step");
    diagnostics_.time = (step_ + 1) * dt_;
    diagnostics_.last_iterations = diagnostics_.blade = 0;
    diagnostics_.block = 0;
    diagnostics_.residual = diagnostics_.scaled_residual = diagnostics_.correction = 0;
    for (std::size_t i = 0; i < layout_.size(); ++i) {
        a0_[i] = (alpha_f_ * acceleration_[i] - alpha_m_ * algorithmic_[i]) * (1 / (1 - alpha_m_));
        predicted_.qd[i] = state_.qd[i] + (dt_ * (1 - gamma_)) * algorithmic_[i] + (gamma_ * dt_) * a0_[i];
        predicted_.q[i] = state_.q[i] + dt_ * state_.qd[i] + (dt_ * dt_ * (.5 - beta_)) * algorithmic_[i] +
                          (beta_ * dt_ * dt_) * a0_[i];
    }
    pending_ = true;
    return predicted_;
}
void GeneralizedAlpha::correct(const AccelerationOperator &op) {
    if (failed_ || !pending_)
        throw std::logic_error("Correction requires a successful prediction");
    try {
        correct_impl(op);
        pending_ = false;
        time_ = (++step_) * dt_;
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void GeneralizedAlpha::correct_impl(const AccelerationOperator &op) {
    state_ = predicted_;
    std::fill(current_.begin(), current_.end(), 0.);
    for (auto &w : work_) {
        w.previous_residual = 0;
        w.age = 0;
    }
    const bool scaled = options_.mode == SolverMode::scaled;
    auto evaluate = [&](std::size_t block, StateView state, std::vector<double> &out) {
        ++diagnostics_.acceleration_evaluations;
        for (std::size_t i = 0; i < state.size; ++i) {
            finite(state.q[i]);
            finite(state.qd[i]);
        }
        std::fill(out.begin(), out.end(), std::numeric_limits<double>::quiet_NaN());
        op.evaluate(block, state, out.data());
        for (double x : out)
            finite(x);
    };
    for (int iteration = 0; iteration < options_.max_iterations; ++iteration) {
        ++diagnostics_.iterations;
        diagnostics_.last_iterations = iteration + 1;
        diagnostics_.max_iterations = std::max(diagnostics_.max_iterations, iteration + 1);
        double largest = 0, worst = 0, worst_scaled = 0;
        for (std::size_t b = 0; b < work_.size(); ++b) {
            diagnostics_.block = b + 1;
            const auto &block = layout_.blocks()[b];
            const auto offset = block.offset, n = block.size;
            auto &w = work_[b];
            evaluate(b, {state_.q.data() + offset, state_.qd.data() + offset, n}, w.f);
            double scaled_residual = 0;
            for (std::size_t i = 0; i < n; ++i) {
                w.rhs[i] = w.f[i] - current_[offset + i];
                finite(w.rhs[i]);
                const double scale = std::max({layout_.dofs()[offset + i].acceleration_scale,
                                               std::abs(w.f[i]), std::abs(current_[offset + i])});
                const double tolerance = options_.residual_absolute + options_.residual_relative * scale;
                if (!positive(tolerance))
                    throw std::runtime_error("Invalid scaled residual tolerance");
                scaled_residual = std::max(scaled_residual, std::abs(w.rhs[i]) / tolerance);
            }
            diagnostics_.residual = vector_norm(w.rhs);
            diagnostics_.scaled_residual = scaled_residual;
            worst = std::max(worst, diagnostics_.residual);
            worst_scaled = std::max(worst_scaled, scaled_residual);
            if (!scaled || !options_.reuse_jacobian || iteration == 0 || w.age >= 2 ||
                scaled_residual >= .9 * w.previous_residual) {
                ++diagnostics_.jacobian_builds;
                w.age = 0;
                for (std::size_t j = 0; j < n; ++j) {
                    const double h = options_.perturbation *
                                     (scaled ? std::max(layout_.dofs()[offset + j].acceleration_scale,
                                                        std::abs(current_[offset + j]))
                                             : 1.);
                    if (!positive(h))
                        throw std::runtime_error("Invalid Jacobian perturbation");
                    std::copy_n(state_.q.data() + offset, n, w.trial.q.begin());
                    std::copy_n(state_.qd.data() + offset, n, w.trial.qd.begin());
                    w.trial.q[j] += beta_prime_ * h;
                    w.trial.qd[j] += gamma_prime_ * h;
                    evaluate(b, {w.trial.q.data(), w.trial.qd.data(), n}, w.fp);
                    for (std::size_t i = 0; i < n; ++i)
                        w.jac[i * n + j] = (i == j ? 1. : 0.) - (w.fp[i] - w.f[i]) / h;
                }
                for (double x : w.jac)
                    finite(x);
            } else
                ++w.age;
            w.previous_residual = scaled_residual;
            w.matrix = w.jac;
            solve(w.matrix, w.rhs, w.delta);
            for (double x : w.delta)
                finite(x);
            const double correction = vector_norm(w.delta);
            finite(correction);
            finite(diagnostics_.residual);
            finite(scaled_residual);
            diagnostics_.correction = correction;
            largest = std::max(largest, correction);
            for (std::size_t i = 0; i < n; ++i) {
                current_[offset + i] += w.delta[i];
                state_.q[offset + i] = predicted_.q[offset + i] + beta_prime_ * current_[offset + i];
                state_.qd[offset + i] = predicted_.qd[offset + i] + gamma_prime_ * current_[offset + i];
                finite(state_.q[offset + i]);
                finite(state_.qd[offset + i]);
            }
        }
        diagnostics_.residual = worst;
        diagnostics_.scaled_residual = worst_scaled;
        diagnostics_.correction = largest;
        if (largest < options_.correction_tolerance && (!scaled || worst_scaled <= 1.))
            break;
        if (iteration + 1 == options_.max_iterations)
            throw std::runtime_error("Structural iteration limit exceeded (reported norms are block maxima)");
    }
    acceleration_ = current_;
    for (std::size_t i = 0; i < layout_.size(); ++i)
        algorithmic_[i] = a0_[i] + ((1 - alpha_f_) / (1 - alpha_m_)) * current_[i];
}
void GeneralizedAlpha::reset() {
    for (auto *v : {&state_.q, &state_.qd, &predicted_.q, &predicted_.qd, &acceleration_, &algorithmic_, &a0_,
                    &current_})
        std::fill(v->begin(), v->end(), 0.);
    time_ = 0;
    step_ = 0;
    pending_ = failed_ = false;
    diagnostics_ = {};
    diagnostics_.mode = options_.mode;
}
} // namespace turbine
