#include "acoustic_quantities.hpp"
#include "aeroacoustics.hpp"
#include <algorithm>
#include <limits>
namespace aeroacoustics {
namespace {
constexpr double pi = 3.14159265358979323846;
bool close(double a, double b) { return std::abs(a - b) <= 1e-12 * std::max(std::abs(a), std::abs(b)); }
double edge(double index) { return 1000. * std::pow(2., index / 3.); }
template <class Q> double total(const BandSpectrum<Q> &s) {
    double result = 0;
    for (double x : s.values())
        result += x;
    if (!std::isfinite(result))
        throw std::overflow_error("Acoustic energy sum overflow");
    return result;
}
template <class Out, class In> Out logarithm(const In &s, double reference) {
    std::vector<double> values;
    for (double x : s.values())
        values.push_back(x == 0 ? -INFINITY : 10. * (std::log10(x) - std::log10(reference)));
    return {s.bands(), std::move(values), s.weighting()};
}
template <class Out, class In> Out linear(const In &s, double reference) {
    std::vector<double> values;
    for (double x : s.values())
        values.push_back(std::pow(10., x / 10. + std::log10(reference)));
    return {s.bands(), std::move(values), s.weighting()};
}
template <class Q> BandSpectrum<Q> weight(const BandSpectrum<Q> &s) {
    if (s.weighting() != Weighting::unweighted)
        throw std::invalid_argument("Spectrum is already A weighted");
    std::vector<double> centers;
    for (const auto &b : s.bands().values())
        centers.push_back(b.center_hz);
    auto correction = a_weighting(centers);
    auto values = s.values();
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] *= std::pow(10., correction[i] / 10.);
    return {s.bands(), std::move(values), Weighting::a};
}
template <class Q> BandSpectrum<Q> merge(const BandSpectrum<Q> &s, FrequencyBands target) {
    const auto &source = s.bands().values();
    std::vector<double> result;
    std::size_t i = 0;
    for (const auto &b : target.values()) {
        if (i == source.size() || !close(source[i].lower_hz, b.lower_hz))
            throw std::invalid_argument("Merge target must cover complete source bands");
        double value = 0, edge_hz = b.lower_hz;
        while (i < source.size() && source[i].upper_hz <= b.upper_hz * (1. + 1e-12)) {
            if (!close(source[i].lower_hz, edge_hz))
                throw std::invalid_argument("Gap inside merged band");
            value += s.values()[i];
            edge_hz = source[i++].upper_hz;
        }
        if (!close(edge_hz, b.upper_hz))
            throw std::invalid_argument("Merge target splits a source band");
        result.push_back(value);
    }
    if (i != source.size())
        throw std::invalid_argument("Merge target omits source bands");
    return {std::move(target), std::move(result), s.weighting()};
}
} // namespace
FrequencyBands::FrequencyBands(std::vector<FrequencyBand> bands) : bands_(std::move(bands)) {
    if (bands_.empty())
        throw std::invalid_argument("Empty frequency bands");
    double last_upper = 0, last_center = 0;
    for (const auto &b : bands_) {
        if (!std::isfinite(b.center_hz) || !std::isfinite(b.lower_hz) || !std::isfinite(b.upper_hz) ||
            b.lower_hz <= 0 || b.upper_hz <= b.lower_hz || b.center_hz <= last_center ||
            b.center_hz < b.lower_hz || b.center_hz > b.upper_hz ||
            (b.lower_hz < last_upper && !close(b.lower_hz, last_upper)))
            throw std::invalid_argument("Frequency bands must be finite, ordered and non-overlapping");
        last_upper = b.upper_hz;
        last_center = b.center_hz;
    }
}
FrequencyBands FrequencyBands::openfast_reference(const std::vector<double> &centers) {
    std::vector<FrequencyBand> bands;
    for (double f : centers) {
        auto it = std::find(openfast_centers_hz.begin(), openfast_centers_hz.end(), f);
        if (it == openfast_centers_hz.end())
            throw std::invalid_argument("Aggregation requires OpenFAST nominal third-octave centers; "
                                        "arbitrary kernel samples are not independent bands");
        const auto k = (it - openfast_centers_hz.begin()) - 20;
        bands.push_back({f, edge(k - .5), edge(k + .5)});
    }
    FrequencyBands result(std::move(bands));
    result.reference_ = true;
    return result;
}
FrequencyBands FrequencyBands::third_octave(int first, std::size_t count) {
    if (count == 0 || count > 10000)
        throw std::invalid_argument("Invalid third-octave band count");
    std::vector<FrequencyBand> bands;
    for (std::size_t i = 0; i < count; ++i) {
        const double k = double(first) + double(i);
        bands.push_back({edge(k), edge(k - .5), edge(k + .5)});
    }
    return FrequencyBands(std::move(bands));
}
bool FrequencyBands::same_as(const FrequencyBands &other) const noexcept {
    if (reference_ != other.reference_ || size() != other.size())
        return false;
    for (std::size_t i = 0; i < size(); ++i) {
        const auto &a = bands_[i], &b = other.bands_[i];
        if (a.center_hz != b.center_hz || a.lower_hz != b.lower_hz || a.upper_hz != b.upper_hz)
            return false;
    }
    return true;
}
ReferenceTnoBandwidth reference_tno_bandwidth(double f) {
    if (!std::isfinite(f) || f <= 0)
        throw std::invalid_argument("Invalid TNO frequency");
    const double omega = 2 * pi * f, ratio = std::pow(2., 1. / 3.);
    const double value = 2 * omega * (std::sqrt(ratio) - 1 / std::sqrt(ratio));
    if (!std::isfinite(value))
        throw std::overflow_error("TNO bandwidth overflow");
    return {value};
}
BandMeanSquarePressure integrate(const PressurePsd &s) {
    auto values = s.values();
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] *= s.bands().values()[i].upper_hz - s.bands().values()[i].lower_hz;
    return {s.bands(), std::move(values), s.weighting()};
}
BandSoundPressureLevel pressure_levels(const BandMeanSquarePressure &s) {
    return logarithm<BandSoundPressureLevel>(s, reference_pressure_pa * reference_pressure_pa);
}
BandMeanSquarePressure mean_square_pressure(const BandSoundPressureLevel &s) {
    return linear<BandMeanSquarePressure>(s, reference_pressure_pa * reference_pressure_pa);
}
BandSoundPowerLevel power_levels(const BandSoundPower &s) {
    return logarithm<BandSoundPowerLevel>(s, reference_sound_power_w);
}
BandSoundPower sound_power(const BandSoundPowerLevel &s) {
    return linear<BandSoundPower>(s, reference_sound_power_w);
}
BandMeanSquarePressure apply_a_weighting(const BandMeanSquarePressure &s) { return weight(s); }
BandSoundPower apply_a_weighting(const BandSoundPower &s) { return weight(s); }
BandMeanSquarePressure merge_bands(const BandMeanSquarePressure &s, FrequencyBands t) {
    return merge(s, std::move(t));
}
BandSoundPower merge_bands(const BandSoundPower &s, FrequencyBands t) { return merge(s, std::move(t)); }
double total_mean_square_pressure(const BandMeanSquarePressure &s) { return total(s); }
double total_sound_power(const BandSoundPower &s) { return total(s); }
namespace {
template <std::size_t... I>
auto typed_mechanisms(const FrequencyBands &bands, Weighting weighting, Mechanisms values,
                      std::index_sequence<I...>) {
    return std::array<BandSoundPressureLevel, mechanism_count>{
        BandSoundPressureLevel(bands, std::move(values[I]), weighting)...};
}
} // namespace
std::array<BandSoundPressureLevel, mechanism_count> band_section_spectrum(const Parameters &p,
                                                                          const Section &s) {
    const auto bands = FrequencyBands::openfast_reference(p.freqlist);
    return typed_mechanisms(bands, p.aweighting ? Weighting::a : Weighting::unweighted,
                            section_spectrum(p, s), std::make_index_sequence<mechanism_count>{});
}
} // namespace aeroacoustics
