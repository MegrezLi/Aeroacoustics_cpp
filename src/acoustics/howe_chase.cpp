#include "source_models.hpp"
#include "trailing_edge.hpp"
#include <algorithm>
#include <cmath>
namespace aeroacoustics {
namespace {
constexpr double pi = 3.14159265358979323846;
void positive(double x, const char *message) {
    if (!std::isfinite(x) || x <= 0)
        throw std::invalid_argument(message);
}
double angular_factor(const Geometry &g) {
    if (!std::isfinite(g.theta) || !std::isfinite(g.phi))
        throw std::invalid_argument("Invalid Howe observer angles");
    positive(g.distance, "Invalid Howe observer distance");
    // Mirror symmetry about the plate; the source paper uses azimuth in [0, pi].
    return std::pow(std::sin(g.theta * pi / 360), 2) * std::abs(std::sin(g.phi * pi / 180));
}
} // namespace
double howe_chase_shape(double q, double h, double wavelength, double delta) {
    positive(delta, "Howe needs positive boundary-layer thickness");
    positive(wavelength, "Howe needs positive wavelength");
    if (!std::isfinite(q) || q < 0 || !std::isfinite(h) || h < 0)
        throw std::invalid_argument("Invalid Howe frequency/half-height");
    // Long-double analytic differentiation of f in Mayer et al., Eqs. 17-18.
    // chi is differentiated at fixed q, h/lambda and h/delta; no finite differences.
    const long double chi = 1.33L, w = q, hh = h, lambda = wavelength, dd = delta;
    const long double w2 = w * w, chi2 = chi * chi, slope = 4 * hh / lambda, a = 1 + slope * slope;
    const long double b = w2 * a + chi2, root = std::sqrt(w2 + chi2);
    long double psi = w2 * a / (b * b);
    if (h > 0 && q > 0) {
        const long double t = lambda / (2 * dd) * root, z = 2 * w * hh / dd;
        // (cosh(t)-cos(z))/sinh(t), evaluated without overflow or small-t subtraction.
        const long double e = std::exp(-t), em = -std::expm1(-t), den = -std::expm1(-2 * t);
        const long double ss = std::pow(std::sin(z / 2), 2);
        const long double H = (em * em + 4 * e * ss) / den;
        const long double Ht = 2 * e * (em * em - 2 * ss * (1 + e * e)) / (den * den);
        const long double c = 64 * (hh / lambda) * (hh / lambda) * (dd / lambda) * w2;
        psi +=
            c / (root * b * b) *
            (H * (1 - chi2 / (2 * root * root) - 2 * chi2 / b) + lambda / (2 * dd) * chi2 / (2 * root) * Ht);
    }
    const double result = static_cast<double>(psi);
    if (!std::isfinite(result) || result < 0)
        throw std::runtime_error("Howe spectrum outside numerical domain");
    return result;
}
double howe_chase_psd(double f, double ue, double delta, double span, const Geometry &g,
                      const Parameters &p) {
    positive(f, "Invalid Howe frequency");
    positive(ue, "Howe needs positive edge speed");
    positive(span, "Invalid Howe span");
    positive(delta, "Invalid Howe boundary layer");
    positive(p.spdsound, "Invalid sound speed");
    positive(p.airdens, "Invalid density");
    if (ue / p.spdsound > .3)
        throw std::invalid_argument("Howe low-Mach model requires Ue/c <= 0.3");
    positive(p.howe.convection_ratio, "Invalid Howe convection ratio");
    positive(p.howe.friction_ratio, "Invalid Howe friction ratio");
    if (p.howe.convection_ratio > 1 || p.howe.friction_ratio >= 1)
        throw std::invalid_argument("Invalid Howe velocity ratios");
    const double uc = p.howe.convection_ratio * ue, ustar = p.howe.friction_ratio * ue;
    const double shape =
        howe_chase_shape(2 * pi * f * delta / uc, p.howe.half_height, p.howe.wavelength, delta);
    // Howe's two-sided angular-frequency PSD: <p^2>=integral_{-inf}^{inf} Phi(omega)domega.
    // Convert once to the one-sided Hz PSD: G(f)=4*pi*Phi(2*pi*f).
    const double psd = 4 * pi * std::pow(p.airdens * ustar, 2) * (span / p.spdsound) *
                       std::pow(delta / g.distance, 2) * (.1553 / pi) * angular_factor(g) * shape;
    if (!std::isfinite(psd) || psd < 0)
        throw std::overflow_error("Invalid Howe PSD");
    return psd;
}
namespace detail {
void prepare_howe(const Parameters &p, const Section &s, HoweSource &out) {
    if (s.speed / p.spdsound > .3 || std::abs(s.alpha_deg) >= std::abs(s.stall_deg))
        throw std::invalid_argument("Howe requires low Mach number and an attached-flow angle below stall");
    if (2 * p.howe.half_height > s.chord)
        throw std::invalid_argument("Howe root-to-tip serration length exceeds section chord");
    const auto bands = FrequencyBands::openfast_reference(p.freqlist);
    for (int side = 0; side < 2; ++side) {
        const double ue = s.speed * s.bl.edge_velocity_ratio[side], delta = s.bl.d99[side];
        if (!std::isfinite(s.bl.cf[side]) || s.bl.cf[side] <= 0)
            throw std::invalid_argument("Howe requires attached turbulent boundary layers with Cf > 0");
        positive(ue, "Howe needs positive input Ue/U on both sides");
        positive(delta, "Howe needs positive d99 on both sides");
        auto &power = out.band_power_at_unit_geometry[side];
        power.resize(bands.size());
        for (std::size_t i = 0; i < bands.size(); ++i) {
            const auto &band = bands.values()[i];
            // At most half a serration phase period per panel, with an explicit work limit.
            const double panels =
                std::max(1., std::ceil((band.upper_hz - band.lower_hz) * 4 * p.howe.half_height /
                                       (p.howe.convection_ratio * ue)));
            if (!std::isfinite(panels) || panels > 1024)
                throw std::invalid_argument("Howe band quadrature too large");
            double value = 0, error = 0;
            for (int k = 0; k < int(panels); ++k) {
                const double lo = band.lower_hz + (band.upper_hz - band.lower_hz) * k / panels;
                const double hi = band.lower_hz + (band.upper_hz - band.lower_hz) * (k + 1) / panels;
                const auto integral = qk61(
                    [&](double f) { return howe_chase_psd(f, ue, delta, s.span, {1, 180, 90}, p); }, lo, hi);
                value += integral.value;
                error += integral.error;
            }
            if (!std::isfinite(value) || value <= 0 || error > 1e-8 * value)
                throw std::runtime_error("Howe band integration failed accuracy check");
            power[i] = value;
        }
    }
}
void emit_howe(const HoweSource &s, const Geometry &g, Spectrum &pressure, Spectrum &suction) {
    const double factor =
        angular_factor(g) / (g.distance * g.distance * reference_pressure_pa * reference_pressure_pa);
    for (int side = 0; side < 2; ++side) {
        auto &out = side == 0 ? suction : pressure;
        const auto &power = s.band_power_at_unit_geometry[side];
        out.resize(power.size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = factor == 0 ? -INFINITY : 10 * std::log10(power[i] * factor);
    }
}
} // namespace detail
} // namespace aeroacoustics
