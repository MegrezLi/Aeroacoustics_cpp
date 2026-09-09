// Lowson/Guidati equations derived from OpenFAST, Apache-2.0; see LICENSE and
// NOTICE.
#include "source_models.hpp"
#include <algorithm>
#include <cmath>

namespace aeroacoustics::detail {
void prepare_inflow(const Parameters &p, const Section &s, double alpha_rad, double ti, InflowSource &out) {
    out.mach = s.speed / p.spdsound;
    const double beta2 = 1. - out.mach * out.mach;
    const double ke = 3. / (4. * p.lturb);
    const double cutoff = (10. * s.speed / pi) / s.chord;
    out.amplitude = std::pow(p.airdens, 2) * std::pow(p.spdsound, 4) * p.lturb * (s.span / 2.);
    out.angle_correction = 10. * std::log10(1. + 9. * std::pow(alpha_rad, 2));
    const auto n = p.freqlist.size();
    out.spectral_power.resize(n);
    out.low_frequency_correction.resize(n);
    out.high_frequency.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double wave = two_pi * p.freqlist[i] / s.speed;
        const double kbar = wave * s.chord / 2., khat = wave / ke;
        out.spectral_power[i] = std::pow(out.mach, 5) * std::pow(ti, 2) * std::pow(khat, 3) *
                                std::pow(1. + std::pow(khat, 2), -7. / 3.);
        const double sears = 1. / ((two_pi * kbar / beta2) + 1. / (1. + 2.4 * kbar / beta2));
        const double lfc = std::max(epsilon, (10. * sears * out.mach) * std::pow(kbar, 2) / beta2);
        out.low_frequency_correction[i] = 10. * log10aa(lfc / (1. + lfc));
        out.high_frequency[i] = p.freqlist[i] > cutoff;
    }
}
void emit_inflow(const InflowSource &s, const Geometry &g, Spectrum &out) {
    const double high = directh_le(s.mach, g.theta, g.phi), low = directl(s.mach, g.theta, g.phi);
    out.resize(s.spectral_power.size());
    if (high <= 0) {
        std::fill(out.begin(), out.end(), 0.);
        return;
    }
    const double amplitude = s.amplitude / std::pow(g.distance, 2);
    for (std::size_t i = 0; i < out.size(); ++i) {
        const double directivity = s.high_frequency[i] ? high : low;
        const double high_level = 10. * log10aa(amplitude * s.spectral_power[i] * directivity) + 78.4;
        out[i] = (high_level + s.angle_correction) + s.low_frequency_correction[i];
    }
}
void prepare_guidati(const Parameters &p, const Section &s, Spectrum &out) {
    const double thickness = s.thickness_1p + s.thickness_10p;
    const double slope = 1.123 * thickness + (5.317 * thickness) * thickness;
    const double frequency_scale = (-slope * two_pi) * s.chord / s.speed;
    const double offset = -slope * 5.;
    out.resize(p.freqlist.size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = frequency_scale * p.freqlist[i] + offset;
}
} // namespace aeroacoustics::detail
