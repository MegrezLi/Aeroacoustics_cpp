// Stable scalar-argument kernel API, used by external callers and Fortran
// comparisons. The implementations share the prepared-source path used by the
// turbine driver.
#include "source_models.hpp"

namespace aeroacoustics {
namespace {
Section section(double alpha, double chord, double speed, double span = 1.) {
    Section s;
    s.alpha_deg = alpha;
    s.chord = chord;
    s.speed = speed;
    s.span = span;
    return s;
}
detail::BpmBoundaryLayer layer(Section &s, const Parameters &p, double delta, double suction, double pressure,
                               double stall) {
    s.stall_deg = stall;
    s.bl.d99[1] = delta;
    s.bl.dstar = {suction, pressure};
    return detail::boundary_layer(p, s);
}
} // namespace
Spectrum lblvs(double alpha, double chord, double speed, double theta, double phi, double span,
               double distance, const Parameters &p, double delta, double suction, double pressure,
               double stall) {
    auto s = section(alpha, chord, speed, span);
    detail::LaminarSource source;
    detail::prepare_laminar(p, s, layer(s, p, delta, suction, pressure, stall), source);
    Spectrum out;
    detail::emit_laminar(source, {distance, theta, phi}, out);
    return out;
}
std::tuple<Spectrum, Spectrum, Spectrum> tblte(double alpha, double chord, double speed, double theta,
                                               double phi, double span, double distance, const Parameters &p,
                                               double delta, double suction, double pressure, double stall) {
    auto s = section(alpha, chord, speed, span);
    detail::TrailingEdgeSource source;
    detail::prepare_trailing_edge(p, s, layer(s, p, delta, suction, pressure, stall), source);
    Spectrum lp, ls, separation;
    detail::emit_trailing_edge(source, {distance, theta, phi}, lp, ls, separation);
    return {std::move(lp), std::move(ls), std::move(separation)};
}
Spectrum blunt(double alpha, double chord, double speed, double theta, double phi, double span,
               double distance, double thickness, double angle, const Parameters &p, double delta,
               double suction, double pressure, double stall) {
    auto s = section(alpha, chord, speed, span);
    s.te_thickness = thickness;
    s.te_angle = angle;
    detail::BluntSource source;
    detail::prepare_blunt(p, s, layer(s, p, delta, suction, pressure, stall), source);
    Spectrum out;
    detail::emit_blunt(source, {distance, theta, phi}, out);
    return out;
}
Spectrum tipnois(double alpha, double lift_ratio, double chord, double speed, double theta, double phi,
                 double distance, const Parameters &p) {
    auto options = p;
    options.alprat = lift_ratio;
    detail::TipSource source;
    detail::prepare_tip(options, section(alpha, chord, speed), source);
    Spectrum out;
    detail::emit_tip(source, {distance, theta, phi}, out);
    return out;
}
Spectrum inflownoise(double alpha_rad, double chord, double speed, double theta, double phi, double span,
                     double distance, double ti, const Parameters &p) {
    detail::InflowSource source;
    detail::prepare_inflow(p, section(0., chord, speed, span), alpha_rad, ti, source);
    Spectrum out;
    detail::emit_inflow(source, {distance, theta, phi}, out);
    return out;
}
Spectrum simple_guidati(double speed, double chord, double thickness10, double thickness1,
                        const Parameters &p) {
    auto s = section(0., chord, speed);
    s.thickness_10p = thickness10;
    s.thickness_1p = thickness1;
    Spectrum out;
    detail::prepare_guidati(p, s, out);
    return out;
}
} // namespace aeroacoustics
