#pragma once
#include "aeroacoustics.hpp"
#include "kernels.hpp"

// Internal source preparation. No observer geometry is stored in these objects.
namespace aeroacoustics::detail {
constexpr double epsilon = 1e-16;
constexpr double pi = 3.14159265358979323846;
constexpr double two_pi = 2 * pi;

struct BpmBoundaryLayer {
    double pressure_thickness, suction_displacement, pressure_displacement;
};
BpmBoundaryLayer boundary_layer(const Parameters &, const Section &);

struct TrailingEdgeSource {
    double mach{}, pressure_displacement{}, suction_displacement{}, span{};
    double k1{}, k2{}, delta_k1{};
    bool separated{};
    Spectrum pressure_shape, suction_shape, separation_shape;
};
void prepare_trailing_edge(const Parameters &, const Section &, const BpmBoundaryLayer &,
                           TrailingEdgeSource &);
void emit_trailing_edge(const TrailingEdgeSource &, const Geometry &, Spectrum &, Spectrum &, Spectrum &);

struct LaminarSource {
    double mach{}, thickness{}, span{};
    Spectrum shape;
};
void prepare_laminar(const Parameters &, const Section &, const BpmBoundaryLayer &, LaminarSource &);
void emit_laminar(const LaminarSource &, const Geometry &, Spectrum &);

struct BluntSource {
    double mach{}, thickness{}, span{}, g4{};
    Spectrum shape, normalization;
};
void prepare_blunt(const Parameters &, const Section &, const BpmBoundaryLayer &, BluntSource &);
void emit_blunt(const BluntSource &, const Geometry &, Spectrum &);

struct TipSource {
    bool zero_alpha{};
    double mach{}, amplitude{};
    Spectrum shape;
};
void prepare_tip(const Parameters &, const Section &, TipSource &);
void emit_tip(const TipSource &, const Geometry &, Spectrum &);

struct InflowSource {
    double mach{}, amplitude{}, angle_correction{};
    Spectrum spectral_power, low_frequency_correction;
    std::vector<bool> high_frequency;
};
// The standalone kernel accepts radians, unlike Section::alpha_deg.
void prepare_inflow(const Parameters &, const Section &, double alpha_rad, double ti, InflowSource &);
void emit_inflow(const InflowSource &, const Geometry &, Spectrum &);
void prepare_guidati(const Parameters &, const Section &, Spectrum &);

struct TnoWorkspace {
    Spectrum wave, height, factor, gauss, decay, exp_gauss, exp_decay, pressure;
    TnoWorkspace();
};
struct TnoSource {
    double mach{}, span{};
    std::array<bool, 2> active{}; // suction, pressure
    std::array<Spectrum, 2> integral;
    Spectrum bandwidth;
};
double integrate_tno(double omega, double lower, double upper, bool suction, double mach,
                     const BoundaryLayer &, const Parameters &, TnoWorkspace &);
void prepare_tno(const Parameters &, const Section &, TnoWorkspace &, TnoSource &);
void emit_tno(const TnoSource &, const Geometry &, Spectrum &pressure, Spectrum &suction);

struct PreparedSection {
    Section section;
    TrailingEdgeSource trailing_edge;
    LaminarSource laminar;
    BluntSource blunt;
    TipSource tip;
    InflowSource inflow;
    TnoSource tno;
    Spectrum guidati;
};
void prepare_section(const Parameters &, const SourceSelection &, const Section &, TnoWorkspace &,
                     PreparedSection &);
void emit_section(const Parameters &, const SourceSelection &, const PreparedSection &,
                  const Geometry &leading, const Geometry &trailing, const Spectrum &weighting, Mechanisms &);
} // namespace aeroacoustics::detail
