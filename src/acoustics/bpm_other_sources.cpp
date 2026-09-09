// BPM equations derived from OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "source_models.hpp"
#include <algorithm>
#include <cmath>

namespace aeroacoustics::detail {
void prepare_laminar(const Parameters &p, const Section &s, const BpmBoundaryLayer &bl, LaminarSource &out) {
    out.mach = s.speed / p.spdsound;
    out.thickness = bl.pressure_thickness;
    out.span = s.span;
    const double re = s.speed * s.chord / p.kinvisc;
    const double st1 = re <= 130000. ? .18 : re <= 400000. ? .001756 * std::pow(re, .3931) : .28;
    const double peak = std::pow(10., -.04 * s.alpha_deg) * st1;
    const double re0 =
        std::pow(10., s.alpha_deg <= 3. ? .215 * s.alpha_deg + 4.978 : .12 * s.alpha_deg + 5.263);
    const double ratio = re / re0;
    double g2;
    if (ratio <= .3237)
        g2 = 77.852 * log10aa(ratio) + 15.328;
    else if (ratio <= .5689)
        g2 = 65.188 * std::log10(ratio) + 9.125;
    else if (ratio <= 1.7579)
        g2 = -114.052 * std::pow(std::log10(ratio), 2);
    else if (ratio <= 3.0889)
        g2 = -65.188 * std::log10(ratio) + 9.125;
    else
        g2 = -77.852 * std::log10(ratio) + 15.328;
    const double g3 = 171.04 - 3.03 * s.alpha_deg;
    out.shape.resize(p.freqlist.size());
    for (std::size_t i = 0; i < p.freqlist.size(); ++i) {
        const double e = (p.freqlist[i] * bl.pressure_thickness / s.speed) / peak;
        double g1;
        if (e <= .5974)
            g1 = 39.8 * log10aa(e) - 11.12;
        else if (e <= .8545)
            g1 = 98.409 * std::log10(e) + 2.;
        else if (e <= 1.17)
            g1 = -5.076 + std::sqrt(2.484 - 506.25 * std::pow(std::log10(e), 2));
        else if (e <= 1.674)
            g1 = -98.409 * std::log10(e) + 2.;
        else
            g1 = -39.8 * std::log10(e) - 11.12;
        out.shape[i] = (g1 + g2) + g3;
    }
}
void emit_laminar(const LaminarSource &s, const Geometry &g, Spectrum &out) {
    const double directivity = directh_te(s.mach, g.theta, g.phi);
    out.resize(s.shape.size());
    if (directivity <= 0) {
        std::fill(out.begin(), out.end(), 0.);
        return;
    }
    const double scale =
        10. * log10aa(s.thickness * std::pow(s.mach, 5) * directivity * s.span / std::pow(g.distance, 2));
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = s.shape[i] + scale;
}

void prepare_blunt(const Parameters &p, const Section &s, const BpmBoundaryLayer &bl, BluntSource &out) {
    out.mach = s.speed / p.spdsound;
    out.thickness = s.te_thickness;
    out.span = s.span;
    const double ratio = s.te_thickness / ((bl.suction_displacement + bl.pressure_displacement) / 2.);
    const double inverse_ratio = 1. / ratio;
    const double peak = ratio >= .2 ? (.212 - .0045 * s.te_angle) /
                                          ((1. + .235 * inverse_ratio) - .0132 * std::pow(inverse_ratio, 2))
                                    : (.1 * ratio + .095) - .00243 * s.te_angle;
    out.g4 = ratio <= 5. ? (17.5 * log10aa(ratio) + 157.5) - 1.114 * s.te_angle : 169.7 - 1.114 * s.te_angle;
    const double modified_ratio = (6.724 * std::pow(ratio, 2) - 4.019 * ratio) + 1.107;
    out.shape.resize(p.freqlist.size());
    out.normalization.resize(p.freqlist.size());
    double cumulative = 0.;
    for (std::size_t i = 0; i < p.freqlist.size(); ++i) {
        const double eta = log10aa((p.freqlist[i] * s.te_thickness / s.speed) / peak);
        const double at14 = g5comp(ratio, eta), at0 = g5comp(modified_ratio, eta);
        out.shape[i] = std::min(at0 + (.0714 * s.te_angle) * (at14 - at0), 0.);
        // Preserve the original prefix sum: each band uses the sum through that
        // band.
        cumulative = std::pow(10., out.shape[i] / 10.) + cumulative;
        const double inverse = cumulative != 0 ? std::max(epsilon, 1. / cumulative) : 1.;
        out.normalization[i] = 10. * std::log10(inverse);
    }
}
void emit_blunt(const BluntSource &s, const Geometry &g, Spectrum &out) {
    const double directivity = directh_te(s.mach, g.theta, g.phi);
    out.resize(s.shape.size());
    if (directivity <= 0) {
        std::fill(out.begin(), out.end(), 0.);
        return;
    }
    const double scale =
        10. * log10aa(std::pow(s.mach, 5.5) * s.thickness * directivity * s.span / std::pow(g.distance, 2));
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = ((s.g4 + s.shape[i]) + scale) - s.normalization[i];
}

void prepare_tip(const Parameters &p, const Section &s, TipSource &out) {
    out.zero_alpha = s.alpha_deg == 0.;
    out.mach = s.speed / p.spdsound;
    out.shape.resize(p.freqlist.size());
    if (out.zero_alpha)
        return;
    const double angle = std::abs(s.alpha_deg) * p.alprat;
    const double length = p.round                 ? .008 * angle * s.chord
                          : std::abs(angle) <= 2. ? (.023 + .0169 * angle) * s.chord
                                                  : (.0378 + .0095 * angle) * s.chord;
    const double maximum_mach = (1. + .036 * angle) * out.mach;
    const double maximum_speed = maximum_mach * p.spdsound;
    out.amplitude = out.mach * out.mach * std::pow(maximum_mach, 3) * std::pow(length, 2);
    for (std::size_t i = 0; i < p.freqlist.size(); ++i) {
        const double st = p.freqlist[i] * length / maximum_speed;
        out.shape[i] = 126. - 30.5 * std::pow(log10aa(st) + .3, 2);
    }
}
void emit_tip(const TipSource &s, const Geometry &g, Spectrum &out) {
    out.resize(s.shape.size());
    if (s.zero_alpha) {
        std::fill(out.begin(), out.end(), 0.);
        return;
    }
    const double term = s.amplitude * directh_te(s.mach, g.theta, g.phi) / std::pow(g.distance, 2);
    const double scale = term != 0 ? 10. * std::log10(term) : 0.;
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = s.shape[i] + scale;
}
} // namespace aeroacoustics::detail
