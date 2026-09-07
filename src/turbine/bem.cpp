#include "turbine/bem.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace turbine {
double prandtl_loss(const BEMOptions &p, const BEMInput &u, double phi) {
    const double s = std::abs(std::sin(phi));
    if (s < 1e-14)
        return 1;
    const auto factor = [&](double c) { return 2 / pi * std::acos(std::min(1., std::exp(-c / s))); };
    return (p.tip_loss ? factor(u.tip_loss_constant) : 1) * (p.hub_loss ? factor(u.hub_loss_constant) : 1);
}
// BEMTUncoupled.inductionFactors0, for BEM_Mod=1 and SkewMomCorr=False.
Induction induction_factors(const BEMOptions &p, const BEMInput &u, double phi, double cn, double ct,
                            double F) {
    Induction y;
    y.phi = phi;
    y.loss = F;
    const double s = std::sin(phi), c = std::cos(phi), solidity = p.blades * u.chord / (2 * pi * u.radius);
    const double k = solidity * cn / (4 * F * s * s);
    y.k = k;
    const bool momentum = (phi > 0 && u.vx >= 0) || (phi < 0 && u.vx < 0);
    if (momentum) {
        if (k <= 2. / 3) {
            y.axial = std::abs(k + 1) < 1e-14 ? -std::copysign(1e6, 1 + k) : k / (1 + k);
            y.valid = k >= -1;
        } else {
            const double t = 2 * F * k, g1 = t - (10. / 9 - F), g2 = t - (4. / 3 - F) * F,
                         g3 = t - (25. / 9 - 2 * F);
            y.axial = std::abs(g3) < 1e-6 ? 1 - .5 / std::sqrt(g2) : (g1 - std::sqrt(std::abs(g2))) / g3;
        }
    } else {
        y.axial = std::abs(k - 1) < 1e-14 ? 1e6 : k / (k - 1);
        y.valid = k > 1;
    }
    if (p.tangential) {
        if (std::abs(c) < 1e-14) {
            y.tangential = -1;
            y.kp = std::copysign(1e6, ct * s) * std::copysign(1., u.vx);
        } else {
            y.kp = solidity * ct / (4 * F * s * c);
            if (u.vx < 0)
                y.kp = -y.kp;
            y.tangential = std::abs(y.kp - 1) < 1e-14 ? std::copysign(1e6, 1 - y.kp) : y.kp / (1 - y.kp);
        }
    }
    if (momentum)
        y.residual = (std::abs(y.axial - 1) < 1e-14 ? 0 : s / (1 - y.axial)) - c * u.vx / u.vy * (1 - y.kp);
    else
        y.residual = s * (1 - k) - c * u.vx / u.vy * (1 - y.kp);
    return y;
}
Induction bem_residual(const BEMOptions &p, const BEMInput &u, const Airfoil &af, double phi) {
    const auto cf = af.at(phi - u.twist);
    const double s = std::sin(phi), c = std::cos(phi);
    const double cn = cf.cl * c + (p.axial_drag ? cf.cd * s : 0),
                 ct = cf.cl * s - (p.tangential_drag ? cf.cd * c : 0);
    return induction_factors(p, u, phi, cn, ct, std::max(.0001, prandtl_loss(p, u, phi)));
}
Induction solve_bem(const BEMOptions &p, const BEMInput &u, const Airfoil &af, double previous_phi) {
    if (!(u.radius > 0 && u.chord > 0))
        throw std::runtime_error("Invalid BEM geometry: radius=" + std::to_string(u.radius) +
                                 ", chord=" + std::to_string(u.chord));
    if ((p.tip_loss && std::abs(u.tip_loss_constant) < 1e-14) ||
        (p.hub_loss && std::abs(u.hub_loss_constant) < 1e-14)) {
        Induction y;
        y.axial = 1;
        y.phi = 0;
        y.loss = 0;
        return y;
    }
    if (std::abs(u.vx) < 1e-12 || std::abs(u.vy) < 1e-12) {
        Induction y;
        y.phi = std::atan2(u.vx, u.vy);
        return y;
    }
    // BEMT.GetSolveRegionOrdering / FindTestRegion: the previous constraint
    // state narrows the bracket. Scanning for the first sign change can select
    // a different root when a stalled polar admits several solutions.
    const double eps = 10 * std::sqrt(std::numeric_limits<double>::epsilon());
    const double direction = u.vx > 0 ? 1. : -1.;
    std::array<std::pair<double, double>, 3> ranges{
        {{eps, pi / 2 - eps}, {-pi / 4, -eps}, {pi / 2 + eps, pi - eps}}};
    if (std::abs(previous_phi) >= pi / 4)
        std::swap(ranges[1], ranges[2]);
    for (auto &range : ranges) {
        range.first *= direction;
        range.second *= direction;
    }
    Induction previous;
    previous.valid = false;
    if (std::abs(previous_phi) > eps && std::abs(std::abs(previous_phi) - pi / 2) > eps) {
        previous = bem_residual(p, u, af, previous_phi);
        if (previous.valid && std::abs(previous.residual) < p.tolerance)
            return previous;
    }
    for (auto range : ranges) {
        double lo = range.first, hi = range.second;
        auto left = bem_residual(p, u, af, lo), right = bem_residual(p, u, af, hi);
        if (left.valid && std::abs(left.residual) < p.tolerance)
            return right.valid && std::abs(right.residual) < std::abs(left.residual) ? right : left;
        if (right.valid && std::abs(right.residual) < p.tolerance)
            return right;
        if (std::signbit(left.residual) == std::signbit(right.residual))
            continue;
        if (previous.valid && previous_phi > std::min(lo, hi) && previous_phi < std::max(lo, hi)) {
            if (std::signbit(left.residual) != std::signbit(previous.residual)) {
                hi = previous_phi;
                right = previous;
            } else {
                lo = previous_phi;
                left = previous;
            }
        }
        // Brent-Dekker in a three-point bracket. The inverse quadratic
        // candidate is evaluated in Lagrange form; interpolation is accepted
        // only when it improves on both the bracket and the previous step.
        auto anchor = left, best = right, opposite = left;
        double step = hi - lo, previous_step = step;
        for (int iteration = 0; iteration < p.max_iterations; ++iteration) {
            if (iteration == 0 || std::signbit(best.residual) == std::signbit(opposite.residual)) {
                opposite = anchor;
                previous_step = best.phi - anchor.phi;
                step = previous_step;
            }
            if (std::abs(opposite.residual) < std::abs(best.residual)) {
                anchor = best;
                best = opposite;
                opposite = anchor;
            }
            const double tolerance = 1e-6 + 2 * std::numeric_limits<double>::epsilon() * std::abs(best.phi);
            const double half = (opposite.phi - best.phi) / 2;
            if ((std::abs(half) <= tolerance || best.residual == 0) && best.valid)
                return best;
            double candidate = (best.phi + opposite.phi) / 2;
            bool interpolate =
                std::abs(previous_step) >= tolerance && std::abs(anchor.residual) > std::abs(best.residual);
            if (interpolate) {
                const double fa = anchor.residual, fb = best.residual, fc = opposite.residual;
                if (anchor.phi == opposite.phi)
                    candidate = best.phi - fb * (best.phi - anchor.phi) / (fb - fa);
                else
                    candidate = anchor.phi * fb * fc / ((fa - fb) * (fa - fc)) +
                                best.phi * fa * fc / ((fb - fa) * (fb - fc)) +
                                opposite.phi * fa * fb / ((fc - fa) * (fc - fb));
                const double proposed = candidate - best.phi, history = previous_step;
                previous_step = step;
                interpolate = std::isfinite(proposed) && proposed * half > 0 &&
                              2 * std::abs(proposed) < 3 * std::abs(half) - tolerance &&
                              std::abs(proposed) < std::abs(history) / 2;
                if (interpolate)
                    step = proposed;
            }
            if (!interpolate) {
                step = half;
                previous_step = half;
            }
            anchor = best;
            best = bem_residual(
                p, u, af, best.phi + (std::abs(step) > tolerance ? step : std::copysign(tolerance, half)));
            if (best.valid && std::abs(best.residual) < p.tolerance)
                return best;
        }
    }
    throw std::runtime_error("BEM could not find a valid inflow-angle root");
}
double skew_axial(double a, double chi0, double ratio, double azimuth, double factor) {
    const double chi = std::remainder((.6 * a + 1) * chi0, 2 * pi);
    const double t = std::abs(chi) > pi / 2 ? std::copysign(1., chi) : std::tan(chi / 2);
    return a * (1 + factor * t * ratio * std::sin(azimuth));
}
} // namespace turbine
