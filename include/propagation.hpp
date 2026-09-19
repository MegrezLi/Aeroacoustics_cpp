#pragma once
#include "aeroacoustics.hpp"
#include <complex>
namespace aeroacoustics {
struct Atmosphere {
    double temperature_k = 293.15, relative_humidity_percent = 70, pressure_pa = 101325;
    void validate() const;
    double absorption_db_per_m(double frequency_hz) const;
};
enum class GroundModel { none, rigid, impedance };
// Vertical, zero-thickness ridge/screen with a horizontal finite top edge.
struct Screen {
    double x1, y1, x2, y2, top_z;
};
struct PropagationOptions {
    bool absorption = true;
    Atmosphere atmosphere;
    GroundModel ground = GroundModel::none;
    double ground_z = 0, sound_speed = 343;
    std::complex<double> normalized_impedance{10, -10};
    std::vector<Screen> screens;
    void validate() const;
};
double knife_edge_loss(double fresnel_parameter);
// Input spectra already include 1/r^2 and source directivity. This layer adds
// only atmospheric loss, an image path, and dominant-screen diffraction.
class OutdoorPropagation {
    PropagationOptions options_;
    FrequencyBands bands_;
    Spectrum absorption_;
    double screen_loss(const Vec3 &, const Vec3 &, double frequency) const;

  public:
    OutdoorPropagation(PropagationOptions, FrequencyBands);
    const PropagationOptions &options() const noexcept { return options_; }
    const FrequencyBands &bands() const noexcept { return bands_; }
    bool has_reflection() const noexcept { return options_.ground != GroundModel::none; }
    Vec3 image_observer(Vec3) const;
    void apply(const Node &, const Vec3 &observer, Mechanisms &direct, const Mechanisms *image) const;
};
} // namespace aeroacoustics
