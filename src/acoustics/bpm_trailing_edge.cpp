// BPM equations derived from OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "source_models.hpp"
#include <algorithm>
#include <cmath>

namespace aeroacoustics::detail {
void prepare_trailing_edge(const Parameters &p, const Section &s, const BpmBoundaryLayer &bl,
                           TrailingEdgeSource &out) {
    const double m = s.speed / p.spdsound;
    const double re = s.speed * s.chord / p.kinvisc;
    const double alpha = s.alpha_deg;
    out.mach = m;
    out.span = s.span;
    out.pressure_displacement = bl.pressure_displacement;
    out.suction_displacement = bl.suction_displacement;
    const double re_pressure = bl.pressure_displacement * s.speed / p.kinvisc;
    const double st1 = 0.02 * std::pow(m, -0.6);
    const double st2 = alpha <= 1.333         ? st1
                       : alpha <= s.stall_deg ? st1 * std::pow(10., 0.0054 * std::pow(alpha - 1.333, 2))
                                              : 4.72 * st1;
    const double st_average = (st1 + st2) / 2.;
    const double min_a = amin(a0comp(re)), max_a = amax(a0comp(re));
    const double min_a3 = amin(a0comp(3. * re)), max_a3 = amax(a0comp(3. * re));
    const double a_ratio = (20. + min_a) / (min_a - max_a);
    const double a3_ratio = (20. + min_a3) / (min_a3 - max_a3);
    const double b0 = re < 95200. ? .3 : re < 857000. ? -4.48e-13 * std::pow(re - 857000., 2) + .56 : .56;
    const double min_b = bmin(b0), max_b = bmax(b0);
    const double b_ratio = (20. + min_b) / (min_b - max_b);
    out.k1 = re < 247000.    ? -4.31 * log10aa(re) + 156.3
             : re <= 800000. ? -9.0 * std::log10(re) + 181.6
                             : 128.5;
    out.delta_k1 = re_pressure <= 5000. ? -alpha * (5.29 - 1.43 * log10aa(re_pressure)) : 0.;
    const double gamma = 27.094 * m + 3.31, beta = 72.65 * m + 10.74;
    const double gamma0 = 23.43 * m + 4.651, beta0 = -34.19 * m - 13.82;
    out.k2 =
        alpha <= gamma0 - gamma ? -1000.
        : alpha <= gamma0 + gamma
            ? std::sqrt(std::pow(beta, 2) - std::pow(beta / gamma, 2) * std::pow(alpha - gamma0, 2)) + beta0
            : -12.;
    out.k2 += out.k1;
    out.separated = alpha >= gamma0 || alpha > s.stall_deg;
    out.pressure_shape.resize(p.freqlist.size());
    out.suction_shape.resize(p.freqlist.size());
    out.separation_shape.resize(p.freqlist.size());
    const auto shape_a = [](double x, double ratio) {
        const double low = amin(x);
        return low + ratio * (amax(x) - low);
    };
    for (std::size_t i = 0; i < p.freqlist.size(); ++i) {
        const double stp = p.freqlist[i] * bl.pressure_displacement / s.speed;
        const double sts = p.freqlist[i] * bl.suction_displacement / s.speed;
        out.pressure_shape[i] = shape_a(log10aa(stp / st1), a_ratio);
        if (!out.separated) {
            out.suction_shape[i] = shape_a(log10aa(sts / st_average), a_ratio);
            const double x = log10aa(sts / st2), low = bmin(x);
            out.separation_shape[i] = low + b_ratio * (bmax(x) - low);
        } else {
            out.suction_shape[i] = 0.;
            out.separation_shape[i] = shape_a(log10aa(sts / st2), a3_ratio);
        }
    }
}

void emit_trailing_edge(const TrailingEdgeSource &s, const Geometry &g, Spectrum &pressure, Spectrum &suction,
                        Spectrum &separation) {
    const double directivity =
        s.separated ? directl(s.mach, g.theta, g.phi) : directh_te(s.mach, g.theta, g.phi);
    const auto scale = [&](double thickness) {
        return 10. *
               log10aa(thickness * std::pow(s.mach, 5) * directivity * s.span / std::pow(g.distance, 2));
    };
    const double pressure_scale = scale(s.pressure_displacement),
                 suction_scale = scale(s.suction_displacement);
    const auto n = s.pressure_shape.size();
    pressure.resize(n);
    suction.resize(n);
    separation.resize(n);
    const auto floor_level = [](double value) { return value < -100. ? -100. : value; };
    for (std::size_t i = 0; i < n; ++i) {
        pressure[i] = floor_level(
            s.separated ? pressure_scale : ((s.pressure_shape[i] + s.k1) - 3. + pressure_scale) + s.delta_k1);
        suction[i] =
            floor_level(s.separated ? suction_scale : (s.suction_shape[i] + s.k1) - 3. + suction_scale);
        separation[i] = floor_level((s.separation_shape[i] + s.k2) + suction_scale);
    }
}
} // namespace aeroacoustics::detail
