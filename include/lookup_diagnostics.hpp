#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace diagnostics {
enum class LookupPolicy { clamp, error };
struct Location {
    double time = 0;
    std::size_t blade = 0, node = 0; // one-based; zero means unspecified
    const char *stage = "unspecified";
};
struct RangeStats {
    std::size_t calls = 0;
    double first_time = 0, last_time = 0, minimum = 0, maximum = 0, lower = 0, upper = 0;
};
using RangeKey = std::tuple<std::string, std::string, std::size_t, std::size_t, std::string>;
struct LookupReport {
    LookupPolicy policy = LookupPolicy::clamp;
    std::map<RangeKey, RangeStats> entries;
    std::size_t calls = 0;
};
// Explicit opt-in, per-thread run context. Nested sessions restore the caller.
inline thread_local LookupReport *active_report = nullptr;
inline thread_local Location active_location{};
class LookupSession {
    LookupReport *previous_;
    Location location_;

  public:
    explicit LookupSession(LookupReport &report) : previous_(active_report), location_(active_location) {
        active_report = &report;
        active_location = {};
    }
    ~LookupSession() {
        active_report = previous_;
        active_location = location_;
    }
    LookupSession(const LookupSession &) = delete;
    LookupSession &operator=(const LookupSession &) = delete;
};
class LookupLocation {
    Location previous_;

  public:
    explicit LookupLocation(Location location) : previous_(active_location) { active_location = location; }
    ~LookupLocation() { active_location = previous_; }
    LookupLocation(const LookupLocation &) = delete;
    LookupLocation &operator=(const LookupLocation &) = delete;
};
inline void check_lookup(const std::string &table, const char *axis, double value, double lower,
                         double upper) {
    if (!std::isfinite(value))
        throw std::invalid_argument("Non-finite lookup input: " + table + " axis=" + axis);
    if (value >= lower && value <= upper)
        return; // exact endpoints are valid
    if (!active_report)
        return; // library callers can install a LookupSession to collect diagnostics
    const auto &where = active_location;
    auto &report = *active_report;
    auto &stats = report.entries[{table, axis, where.blade, where.node, where.stage}];
    if (!stats.calls) {
        stats.first_time = stats.last_time = where.time;
        stats.minimum = stats.maximum = value;
        stats.lower = lower;
        stats.upper = upper;
    }
    ++stats.calls;
    ++report.calls;
    stats.first_time = std::min(stats.first_time, where.time);
    stats.last_time = std::max(stats.last_time, where.time);
    stats.minimum = std::min(stats.minimum, value);
    stats.maximum = std::max(stats.maximum, value);
    if (report.policy == LookupPolicy::error) {
        std::ostringstream message;
        message.precision(17);
        message << "Lookup out of range: table=" << table << " axis=" << axis << " value=" << value
                << " range=[" << lower << ',' << upper << "] time=" << where.time << " blade=" << where.blade
                << " node=" << where.node << " stage=" << where.stage;
        throw std::out_of_range(message.str());
    }
}
} // namespace diagnostics
