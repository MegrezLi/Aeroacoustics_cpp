// BPM boundary-layer correlations from OpenFAST, Apache-2.0; see NOTICE.
#include "source_models.hpp"
#include <algorithm>
#include <cmath>

namespace aeroacoustics {
double log10aa(double x) { return std::log10(std::max(detail::epsilon, x)); }

std::tuple<double, double, double> thick(double chord, double re, double alpha, const Parameters &p,
                                         double stall) {
    const double log_re = log10aa(re);
    double delta0 = std::pow(10., (1.6569 - .9045 * log_re) + .0596 * std::pow(log_re, 2)) * chord;
    if (p.itrip != 0)
        delta0 = std::pow(10., (1.892 - .9045 * log_re) + .0596 * std::pow(log_re, 2)) * chord;
    if (p.itrip == 2)
        delta0 *= .6;
    const double pressure_thickness = std::pow(10., -.04175 * alpha + .00106 * std::pow(alpha, 2)) * delta0;
    double displacement0;
    if (p.itrip != 0) {
        displacement0 = re <= 300000.
                            ? (.0601 * std::pow(re, -.114)) * chord
                            : std::pow(10., (3.411 - 1.5397 * log_re) + .1059 * std::pow(log_re, 2)) * chord;
        if (p.itrip == 2)
            displacement0 *= .6;
    } else {
        displacement0 = std::pow(10., (3.0187 - 1.5397 * log_re) + .1059 * std::pow(log_re, 2)) * chord;
    }
    const double pressure_displacement =
        std::pow(10., -.0432 * alpha + .00113 * std::pow(alpha, 2)) * displacement0;
    double suction_displacement;
    if (p.itrip == 1) {
        suction_displacement = alpha <= 5.      ? std::pow(10., .0679 * alpha) * displacement0
                               : alpha <= stall ? (.381 * std::pow(10., .1516 * alpha)) * displacement0
                                                : (14.296 * std::pow(10., .0258 * alpha)) * displacement0;
    } else {
        suction_displacement = alpha <= 7.5     ? std::pow(10., .0679 * alpha) * displacement0
                               : alpha <= stall ? (.0162 * std::pow(10., .3066 * alpha)) * displacement0
                                                : (52.42 * std::pow(10., .0258 * alpha)) * displacement0;
    }
    return {pressure_thickness, suction_displacement, pressure_displacement};
}

namespace detail {
BpmBoundaryLayer boundary_layer(const Parameters &p, const Section &s) {
    if (p.x_blmethod == 2)
        return {s.bl.d99[1], s.bl.dstar[0], s.bl.dstar[1]};
    const auto [thickness, suction, pressure] =
        thick(s.chord, s.speed * s.chord / p.kinvisc, s.alpha_deg, p, s.stall_deg);
    return {thickness, suction, pressure};
}
} // namespace detail
} // namespace aeroacoustics
