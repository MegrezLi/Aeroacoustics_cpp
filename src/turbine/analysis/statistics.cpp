#include "turbine/metrics.hpp"
#include <algorithm>
#include <complex>
#include <numeric>
namespace turbine {
namespace {
double level(double p) { return p == 0 ? -INFINITY : 10 * std::log10(p); }
} // namespace
LevelStatistics level_statistics(const std::vector<double> &t, const std::vector<double> &p) {
    if (t.size() < 2 || t.size() != p.size())
        throw std::invalid_argument("Statistics require two or more samples");
    double energy = 0;
    for (std::size_t i = 0; i < t.size(); ++i) {
        if (!std::isfinite(t[i]) || !std::isfinite(p[i]) || p[i] < 0 || (i && t[i] <= t[i - 1]))
            throw std::invalid_argument("Invalid statistical time/energy samples");
        if (i)
            energy += (t[i] - t[i - 1]) * (.5 * p[i] + .5 * p[i - 1]);
    }
    const double duration = t.back() - t.front();
    auto quantile = [&](double fraction) {
        double lo = *std::min_element(p.begin(), p.end()), hi = *std::max_element(p.begin(), p.end());
        for (int iteration = 0; iteration < 80 && lo < hi; ++iteration) {
            double q = .5 * lo + .5 * hi, cumulative = 0;
            for (std::size_t i = 1; i < t.size(); ++i) {
                const double a = std::min(p[i - 1], p[i]), b = std::max(p[i - 1], p[i]);
                cumulative +=
                    (t[i] - t[i - 1]) * (a == b ? double(q >= a) : std::clamp((q - a) / (b - a), 0., 1.));
            }
            if (cumulative < fraction * duration)
                lo = q;
            else
                hi = q;
        }
        return level(.5 * lo + .5 * hi);
    };
    if (!std::isfinite(energy))
        throw std::overflow_error("Statistical energy overflow");
    return {duration, level(energy / duration), quantile(.95), quantile(.5), quantile(.05)};
}
Modulation modulation(const std::vector<double> &levels, double dt, double low, double high) {
    if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(low) || !std::isfinite(high) || low <= 0 ||
        high < low || high >= .5 / dt)
        throw std::invalid_argument("AM frequency range must lie below sampling Nyquist");
    Modulation result;
    if (levels.size() < 8)
        return result;
    for (double x : levels)
        if (!std::isfinite(x))
            return result; // silence is unresolved, not zero dB
    const std::size_t n = levels.size();
    const double duration = n * dt, mean = std::accumulate(levels.begin(), levels.end(), 0.) / n;
    auto coefficient = [&](std::size_t k) {
        std::complex<double> value{};
        for (std::size_t j = 0; j < n; ++j)
            value += (levels[j] - mean) * std::polar(1., -2 * pi * k * j / n);
        return value / double(n);
    };
    const auto first = std::max<std::size_t>(1, std::ceil(low * duration)),
               last = std::size_t(std::floor(high * duration));
    if (last < first)
        return result;
    std::size_t peak = first;
    double best = 0;
    for (std::size_t k = first; k <= last; ++k) {
        double a = std::norm(coefficient(k));
        if (a > best) {
            best = a;
            peak = k;
        }
    }
    if (best < 1e-20) {
        result.resolved = true;
        return result;
    }
    std::vector<double> reconstruction(n, 0.);
    for (std::size_t h = 1; h <= 3 && h * peak < n / 2; ++h) {
        const auto c = coefficient(h * peak);
        for (std::size_t j = 0; j < n; ++j)
            reconstruction[j] += 2 * std::real(c * std::polar(1., 2 * pi * h * peak * j / n));
    }
    result.resolved = true;
    result.frequency_hz = peak / duration;
    result.harmonic_depth_db = *std::max_element(reconstruction.begin(), reconstruction.end()) -
                               *std::min_element(reconstruction.begin(), reconstruction.end());
    std::sort(reconstruction.begin(), reconstruction.end());
    auto percentile = [&](double f) {
        double q = f * (n - 1);
        auto a = std::size_t(q);
        return reconstruction[a] + (q - a) * (reconstruction[std::min(a + 1, n - 1)] - reconstruction[a]);
    };
    result.percentile_depth_db = percentile(.95) - percentile(.05);
    return result;
}
double apparent_sound_power(double lp, double r, double correction) {
    if (std::isnan(lp) || lp == INFINITY || !std::isfinite(r) || r <= 0 || !std::isfinite(correction))
        throw std::invalid_argument("Invalid apparent sound power geometry/level");
    return lp + 10 * std::log10(4 * pi) + 20 * std::log10(r) - correction;
}
void ArrivalSeries::append(ArrivalSample s) {
    if (!std::isfinite(s.reception) || !std::isfinite(s.emission) || !std::isfinite(s.wind) || s.wind < 0 ||
        s.power.empty() || !std::isfinite(s.frequency_hz) || s.frequency_hz < 0 ||
        !std::isfinite(s.distance_m) || s.distance_m < 0)
        throw std::invalid_argument("Invalid arrival sample");
    for (double p : s.power)
        if (!std::isfinite(p) || p < 0)
            throw std::invalid_argument("Invalid arrival energy");
    if (!samples_.empty() &&
        (s.reception <= samples_.back().reception || s.emission <= samples_.back().emission ||
         s.power.size() != samples_.back().power.size()))
        throw std::invalid_argument("Nonmonotone/subsonic arrival or inconsistent spectrum");
    samples_.push_back(std::move(s));
}
namespace {
std::size_t interval(const std::vector<ArrivalSample> &s, double t) {
    if (s.size() < 2 || !std::isfinite(t) || t < s.front().reception || t > s.back().reception)
        throw std::out_of_range("Receiver time outside available emission history");
    auto it = std::upper_bound(s.begin(), s.end(), t,
                               [](double x, const ArrivalSample &p) { return x < p.reception; });
    return std::min<std::size_t>(std::size_t(it - s.begin()), s.size() - 1) - 1;
}
} // namespace
ArrivalSample ArrivalSeries::at(double t) const {
    auto i = interval(samples_, t);
    const auto &a = samples_[i], &b = samples_[i + 1];
    double f = (t - a.reception) / (b.reception - a.reception);
    ArrivalSample r{t, a.emission + f * (b.emission - a.emission), a.wind + f * (b.wind - a.wind), a.power};
    r.frequency_hz = (1 - f) * a.frequency_hz + f * b.frequency_hz;
    r.distance_m = (1 - f) * a.distance_m + f * b.distance_m;
    for (std::size_t j = 0; j < r.power.size(); ++j)
        r.power[j] = (1 - f) * a.power[j] + f * b.power[j];
    return r;
}
double ArrivalSeries::doppler(double t) const {
    auto i = interval(samples_, t);
    return (samples_[i + 1].emission - samples_[i].emission) /
           (samples_[i + 1].reception - samples_[i].reception);
}
} // namespace turbine
