#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace aeroacoustics {
inline constexpr double reference_pressure_pa = 2e-5;
inline constexpr double reference_sound_power_w = 1e-12;
inline constexpr std::array<double, 34> openfast_centers_hz{
    10,  12.5, 16,  20,   25,   31.5, 40,   50,   63,   80,   100,  125,  160,  200,   250,   315,   400,
    500, 630,  800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000};
enum class Weighting { unweighted, a };
struct FrequencyBand {
    double center_hz, lower_hz, upper_hz;
};
class FrequencyBands {
    std::vector<FrequencyBand> bands_;
    bool reference_ = false;

  public:
    explicit FrequencyBands(std::vector<FrequencyBand>);
    // Binary third-octave edges; rounded OpenFAST centers remain evaluation/label frequencies.
    static FrequencyBands openfast_reference(const std::vector<double> &centers);
    static FrequencyBands third_octave(int first_index, std::size_t count);
    const std::vector<FrequencyBand> &values() const noexcept { return bands_; }
    std::size_t size() const noexcept { return bands_.size(); }
    bool is_openfast_reference() const noexcept { return reference_; }
    bool same_as(const FrequencyBands &) const noexcept;
};
// This is the legacy TNO multiplier, not a bandwidth in Hz or an SI PSD integral.
struct ReferenceTnoBandwidth {
    double value;
};
ReferenceTnoBandwidth reference_tno_bandwidth(double center_hz);

struct PressurePsdTag {
    static constexpr bool logarithmic = false;
    static constexpr const char *unit = "Pa^2/Hz";
};
struct MeanSquarePressureTag {
    static constexpr bool logarithmic = false;
    static constexpr const char *unit = "Pa^2";
};
struct SoundPressureLevelTag {
    static constexpr bool logarithmic = true;
    static constexpr const char *unit = "dB re 20 uPa";
};
struct SoundPowerTag {
    static constexpr bool logarithmic = false;
    static constexpr const char *unit = "W";
};
struct SoundPowerLevelTag {
    static constexpr bool logarithmic = true;
    static constexpr const char *unit = "dB re 1 pW";
};
template <class Quantity> class BandSpectrum {
    FrequencyBands bands_;
    std::vector<double> values_;
    Weighting weighting_;

  public:
    BandSpectrum(FrequencyBands bands, std::vector<double> values,
                 Weighting weighting = Weighting::unweighted)
        : bands_(std::move(bands)), values_(std::move(values)), weighting_(weighting) {
        if (weighting != Weighting::unweighted && weighting != Weighting::a)
            throw std::invalid_argument("Unknown acoustic weighting");
        if (bands_.size() != values_.size())
            throw std::invalid_argument("Band spectrum shape mismatch");
        for (double x : values_)
            if constexpr (Quantity::logarithmic) {
                if (std::isnan(x) || x == INFINITY)
                    throw std::invalid_argument("Invalid acoustic level");
            } else if (!std::isfinite(x) || x < 0)
                throw std::invalid_argument("Invalid acoustic energy/density");
    }
    const FrequencyBands &bands() const noexcept { return bands_; }
    const std::vector<double> &values() const noexcept { return values_; }
    Weighting weighting() const noexcept { return weighting_; }
    static constexpr const char *unit() noexcept { return Quantity::unit; }
};
using PressurePsd = BandSpectrum<PressurePsdTag>;
using BandMeanSquarePressure = BandSpectrum<MeanSquarePressureTag>;
using BandSoundPressureLevel = BandSpectrum<SoundPressureLevelTag>;
using BandSoundPower = BandSpectrum<SoundPowerTag>;
using BandSoundPowerLevel = BandSpectrum<SoundPowerLevelTag>;
BandMeanSquarePressure integrate(const PressurePsd &);
BandSoundPressureLevel pressure_levels(const BandMeanSquarePressure &);
BandMeanSquarePressure mean_square_pressure(const BandSoundPressureLevel &);
BandSoundPowerLevel power_levels(const BandSoundPower &);
BandSoundPower sound_power(const BandSoundPowerLevel &);
BandMeanSquarePressure apply_a_weighting(const BandMeanSquarePressure &);
BandSoundPower apply_a_weighting(const BandSoundPower &);
// Merge whole contiguous bands only; no interpolation or invented sub-band shape.
BandMeanSquarePressure merge_bands(const BandMeanSquarePressure &, FrequencyBands target);
BandSoundPower merge_bands(const BandSoundPower &, FrequencyBands target);
double total_mean_square_pressure(const BandMeanSquarePressure &); // Pa^2
double total_sound_power(const BandSoundPower &);                  // W
struct AcousticMetadata {
    FrequencyBands bands;
    Weighting weighting;
    // Legacy result.power is <p^2>/(20 uPa)^2, never watts or PSD.
    const char *linear_quantity() const noexcept { return "relative_mean_square_pressure"; }
    const char *level_unit() const noexcept { return weighting == Weighting::a ? "dBA" : "dB"; }
};
} // namespace aeroacoustics
