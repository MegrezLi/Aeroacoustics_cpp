#include "turbine/surface.hpp"
#include <algorithm>
namespace turbine {
SurfaceData::SurfaceData(const std::filesystem::path &path)
    : airfoil(0),
      boundary_layer(aeroacoustics::BLTable::read(InputFile(path).file("BoundaryLayer").string())),
      prepared(boundary_layer) {
    InputFile f(path);
    const int id = f.integer("AirfoilID");
    if (id < 1)
        throw std::invalid_argument("Surface AirfoilID must be one based");
    airfoil = std::size_t(id - 1);
    name = f.value("State");
    provenance = f.value("Provenance");
    uncertainty_note = f.value("UncertaintyNote");
    alpha_min = f.number("AlphaMinDeg");
    alpha_max = f.number("AlphaMaxDeg");
    re_min = f.number("ReMin");
    re_max = f.number("ReMax");
    transition_suction = f.number("TransitionSuction");
    transition_pressure = f.number("TransitionPressure");
    roughness_m = f.number("RoughnessM");
    erosion_m = f.number("ErosionM");
    relative_uncertainty = f.number("RelativeUncertainty");
    te_thickness_m = f.number("TEThicknessM");
    te_angle_deg = f.number("TEAngleDeg");
    if (name.empty() || provenance.empty() || uncertainty_note.empty() || alpha_min >= alpha_max ||
        re_min <= 0 || re_min >= re_max || alpha_min < boundary_layer.aoa.front() ||
        alpha_max > boundary_layer.aoa.back() || re_min < boundary_layer.reynolds.front() ||
        re_max > boundary_layer.reynolds.back() || transition_suction < 0 || transition_suction > 1 ||
        transition_pressure < 0 || transition_pressure > 1 || roughness_m < 0 || erosion_m < 0 ||
        relative_uncertainty < 0 || relative_uncertainty > 1 || te_thickness_m <= 0 || te_angle_deg < 0 ||
        te_angle_deg > 14)
        throw std::invalid_argument("Invalid surface metadata/coverage: " + path.string());
    for (const auto &v : boundary_layer.values)
        for (int side = 0; side < 2; ++side)
            if (v[side] <= 0 || v[2 + side] <= 0 || v[4 + side] < v[2 + side] || v[6 + side] < 0)
                throw std::invalid_argument("Surface BL requires Ue>0, 0<dstar<=d99, Cf>=0");
    if (f.value("Polar") != "none") {
        polar.emplace(f.file("Polar"));
        if (alpha_min < polar->alpha.front() / deg || alpha_max > polar->alpha.back() / deg)
            throw std::invalid_argument("Surface polar does not cover declared AoA interval");
        for (const auto &v : polar->coefficients)
            if (!std::isfinite(v.cl) || !std::isfinite(v.cd) || !std::isfinite(v.cm) || v.cd < 0)
                throw std::invalid_argument("Invalid surface polar coefficients");
    }
}
aeroacoustics::BoundaryLayer SurfaceData::at(double alpha, double re, double chord) const {
    if (!std::isfinite(alpha) || !std::isfinite(re) || alpha < alpha_min || alpha > alpha_max ||
        re < re_min || re > re_max)
        throw std::out_of_range("Surface dataset " + name + " outside validated AoA/Re range: alpha=" +
                                std::to_string(alpha) + ", Re=" + std::to_string(re));
    return prepared.interpolate(alpha, re, chord);
}
std::shared_ptr<const SurfaceSet> SurfaceSet::read(const std::filesystem::path &path) {
    InputFile f(path);
    int n = f.integer("NumStates");
    if (n < 1 || n > 10000)
        throw std::invalid_argument("Invalid surface dataset count");
    auto result = std::make_shared<SurfaceSet>();
    for (const auto &file : f.files_after("StateFiles", std::size_t(n))) {
        SurfaceData data(file);
        if (result->find(data.airfoil))
            throw std::invalid_argument("Duplicate surface AirfoilID");
        result->data_.push_back(std::move(data));
    }
    return result;
}
const SurfaceData *SurfaceSet::find(std::size_t id) const noexcept {
    for (const auto &d : data_)
        if (d.airfoil == id)
            return &d;
    return nullptr;
}
Case SurfaceSet::apply(const Case &original) const {
    Case c = original;
    for (const auto &d : data_) {
        if (d.airfoil >= c.airfoils.size())
            throw std::invalid_argument("Surface AirfoilID outside case");
        if (d.polar)
            c.airfoils[d.airfoil] = *d.polar;
    }
    c.validate_scope();
    return c;
}
} // namespace turbine
