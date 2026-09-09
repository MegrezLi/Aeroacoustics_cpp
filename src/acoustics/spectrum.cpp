// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include "source_models.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aeroacoustics {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double silence = -std::numeric_limits<double>::infinity();
} // namespace
void validate(const Parameters &p) {
    if (p.freqlist.empty())
        throw std::invalid_argument("Empty frequency list");
    double last = 0;
    for (double f : p.freqlist) {
        if (!std::isfinite(f) || f <= last)
            throw std::invalid_argument("Frequencies must be finite, positive and increasing");
        last = f;
    }
    for (double v : {p.spdsound, p.kinvisc, p.airdens, p.lturb, p.alprat})
        if (!std::isfinite(v) || v <= 0)
            throw std::invalid_argument("Positive physical parameters required");
    for (double v : {p.ti, p.avgv})
        if (!std::isfinite(v) || v < 0)
            throw std::invalid_argument("Nonnegative TI and mean speed required");
    if (p.timod < 0 || p.timod > 2 || p.tbltemod < 0 || p.tbltemod > 2 || p.x_blmethod < 1 ||
        p.x_blmethod > 2 || p.itrip < 0 || p.itrip > 2 || p.lammod < 0 || p.lammod > 1 || p.tipmod < 0 ||
        p.tipmod > 1 || p.bluntmod < 0 || p.bluntmod > 1)
        throw std::invalid_argument("Unsupported acoustic switch");
}
Spectrum a_weighting(const Spectrum &frequencies) {
    Spectrum result;
    result.reserve(frequencies.size());
    for (double f : frequencies) {
        const double f2 = f * f, f4 = f2 * f2;
        result.push_back(
            10 * std::log10(1.562339 * f4 / ((f2 + std::pow(107.65265, 2)) * (f2 + std::pow(737.86223, 2)))) +
            10 * std::log10(
                     2.242881e16 * f4 /
                     (std::pow(f2 + std::pow(20.598997, 2), 2) * std::pow(f2 + std::pow(12194.22, 2), 2))));
    }
    return result;
}
double db_sum(const Spectrum &levels) {
    if (levels.empty())
        return silence;
    const double maximum = *std::max_element(levels.begin(), levels.end());
    if (maximum == silence)
        return silence;
    double total = 0.;
    for (double level : levels)
        total += std::pow(10., (level - maximum) / 10.);
    return maximum + 10 * std::log10(total);
}
namespace detail {
void prepare_section(const Parameters &p, const Section &input, TnoWorkspace &work, PreparedSection &out) {
    for (double value : {input.speed, input.chord, input.span, input.alpha_deg, input.stall_deg})
        if (!std::isfinite(value))
            throw std::invalid_argument("Section inputs must be finite");
    if (input.speed < 0 || input.speed >= p.spdsound || input.chord <= 0 || input.span <= 0)
        throw std::invalid_argument("Invalid section speed or geometry");
    if (p.x_blmethod == 2 &&
        (input.bl.dstar[0] <= 0 || input.bl.dstar[1] <= 0 || input.bl.d99[0] <= 0 || input.bl.d99[1] <= 0))
        throw std::invalid_argument("Tabulated BL/TNO needs dimensional boundary-layer data");
    out.section = input;
    auto &s = out.section;
    s.speed = std::max(input.speed, .1);
    s.alpha_deg = input.alpha_deg - 360 * std::floor((input.alpha_deg + 180) / 360);
    BpmBoundaryLayer bl{};
    if ((p.lammod && p.itrip == 0) || p.tbltemod || p.bluntmod)
        bl = boundary_layer(p, s);
    if (p.lammod && p.itrip == 0)
        prepare_laminar(p, s, bl, out.laminar);
    if (p.tbltemod)
        prepare_trailing_edge(p, s, bl, out.trailing_edge);
    if (p.tbltemod == 2) {
        // Keep the reference section driver convention; the direct TNO kernel
        // still accepts caller-supplied edge velocity ratios.
        auto tno_section = s;
        tno_section.bl.edge_velocity_ratio = {1., 1.};
        prepare_tno(p, tno_section, work, out.tno);
    }
    if (p.bluntmod) {
        if (s.te_thickness <= 0 || s.te_angle < 0 || s.te_angle > 14)
            throw std::invalid_argument("Invalid bluntness geometry");
        prepare_blunt(p, s, bl, out.blunt);
    }
    if (p.tipmod && s.is_tip)
        prepare_tip(p, s, out.tip);
    if (p.timod) {
        const double ti = s.ti_section < 0 ? p.ti * p.avgv / s.speed : s.ti_section;
        if (!std::isfinite(ti) || ti < 0)
            throw std::invalid_argument("Invalid inflow noise input");
        prepare_inflow(p, s, s.alpha_deg * detail::pi / 180, ti, out.inflow);
        if (p.timod == 2)
            prepare_guidati(p, s, out.guidati);
    }
}
void emit_section(const Parameters &p, const PreparedSection &s, const Geometry &leading,
                  const Geometry &trailing, const Spectrum &weighting, Mechanisms &out) {
    for (double value : {trailing.distance, trailing.theta, trailing.phi})
        if (!std::isfinite(value))
            throw std::invalid_argument("Section geometry must be finite");
    if (trailing.distance <= 0)
        throw std::invalid_argument("Invalid observer distance");
    for (auto &v : out) {
        v.resize(p.freqlist.size());
        std::fill(v.begin(), v.end(), silence);
    }
    if (p.lammod && p.itrip == 0)
        emit_laminar(s.laminar, trailing, out[0]);
    if (p.tbltemod)
        emit_trailing_edge(s.trailing_edge, trailing, out[1], out[2], out[3]);
    if (p.tbltemod == 2)
        emit_tno(s.tno, trailing, out[1], out[2]);
    if (p.bluntmod)
        emit_blunt(s.blunt, trailing, out[4]);
    if (p.tipmod && s.section.is_tip)
        emit_tip(s.tip, trailing, out[5]);
    if (p.timod) {
        if (!std::isfinite(leading.distance) || !std::isfinite(leading.theta) ||
            !std::isfinite(leading.phi) || leading.distance <= 0)
            throw std::invalid_argument("Invalid leading-edge geometry");
        emit_inflow(s.inflow, leading, out[6]);
        if (p.timod == 2)
            for (std::size_t i = 0; i < p.freqlist.size(); ++i)
                out[6][i] += s.guidati[i] + 10.;
    }
    if (p.aweighting)
        for (auto &v : out)
            for (std::size_t i = 0; i < v.size(); ++i)
                v[i] += weighting[i];
}
} // namespace detail

Mechanisms section_spectrum(const Parameters &input, const Section &s) {
    validate(input);
    auto parameters = input;
    if (parameters.tbltemod == 2)
        parameters.x_blmethod = 2;
    detail::TnoWorkspace work;
    detail::PreparedSection source;
    detail::prepare_section(parameters, s, work, source);
    const auto weighting = parameters.aweighting ? a_weighting(parameters.freqlist) : Spectrum{};
    Mechanisms out;
    detail::emit_section(parameters, source, s.leading.value_or(s.trailing), s.trailing, weighting, out);
    return out;
}
} // namespace aeroacoustics
