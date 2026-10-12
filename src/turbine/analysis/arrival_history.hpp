#pragma once
#include "turbine/metrics.hpp"
#include <algorithm>

namespace turbine::detail {
// Internal broadband storage. Public ArrivalSeries keeps its existing owning API.
// Two contiguous arrays avoid one allocation per source sample. Copies own their data.
class ArrivalHistory {
  public:
    struct Sample { double reception, emission, wind; };

  private:
    std::vector<Sample> samples_;
    std::vector<double> power_;
    std::size_t frequencies_ = 0;

  public:
    const std::vector<Sample> &samples() const noexcept { return samples_; }
    double power(std::size_t sample, std::size_t frequency) const {
        return power_[sample * frequencies_ + frequency];
    }
    void append(double reception, double emission, double wind, const std::vector<double> &power) {
        if (!std::isfinite(reception) || !std::isfinite(emission) || !std::isfinite(wind) || wind < 0 ||
            power.empty() || (!samples_.empty() && (power.size() != frequencies_ ||
            reception <= samples_.back().reception || emission <= samples_.back().emission)))
            throw std::invalid_argument("Invalid broadband arrival sample");
        for (double value : power)
            if (!std::isfinite(value) || value < 0)
                throw std::invalid_argument("Invalid arrival energy");
        const auto previous = power_.size();
        power_.insert(power_.end(), power.begin(), power.end());
        try {
            samples_.push_back({reception, emission, wind});
        } catch (...) {
            power_.resize(previous);
            throw;
        }
        frequencies_ = power.size();
    }
    // Cursor belongs to one traversal, not to the history, preserving const/thread safety.
    void add_power(double time, std::vector<double> &out, std::size_t &cursor) const {
        if (out.size() != frequencies_ || samples_.size() < 2 || !std::isfinite(time) ||
            time < samples_.front().reception || time > samples_.back().reception)
            throw std::out_of_range("Invalid broadband interpolation request");
        if (cursor >= samples_.size() - 1 || time < samples_[cursor].reception) {
            const auto it = std::upper_bound(samples_.begin(), samples_.end(), time,
                [](double t, const Sample &s) { return t < s.reception; });
            cursor = std::min<std::size_t>(it - samples_.begin(), samples_.size() - 1) - 1;
        } else {
            while (cursor + 2 < samples_.size() && time >= samples_[cursor + 1].reception)
                ++cursor;
        }
        const auto &a = samples_[cursor], &b = samples_[cursor + 1];
        const double f = (time - a.reception) / (b.reception - a.reception);
        for (std::size_t j = 0; j < frequencies_; ++j)
            out[j] += (1 - f) * power(cursor, j) + f * power(cursor + 1, j);
    }
};
} // namespace turbine::detail
