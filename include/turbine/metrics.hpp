#pragma once
#include "propagation.hpp"
#include "turbine/input.hpp"
namespace turbine {
struct ToneSource {
    std::string name, provenance;
    std::size_t blade = 0, node = 0; // zero based, attaches to aerodynamic center
    double frequency_hz = 0, rotor_order = 0, level_db = 0, reference_distance = 1;
};
struct MetricsOptions {
    double start = 0, end = 20, receiver_dt = .1, am_window = 10, am_min_hz = .2, am_max_hz = 1;
    double wind_bin_width = 1;
    bool retarded_time = true, apparent_power = false;
    // Optional background is added to the predicted turbine signal. It is never subtracted twice.
    std::optional<double> background_laeq_db;
    std::vector<Vec3> observers; // optional replacement observer map; one-based IDs in output
    std::vector<ToneSource> tones;
    std::size_t max_values = 32000000;
    static MetricsOptions read(const std::filesystem::path &);
    void validate() const;
};
struct LevelStatistics {
    double duration = 0, leq_db = -INFINITY, l5_db = -INFINITY, l50_db = -INFINITY, l95_db = -INFINITY;
};
// Linear interpolation of mean-square pressure, exact trapezoidal time integral.
LevelStatistics level_statistics(const std::vector<double> &time, const std::vector<double> &power);
struct Modulation {
    bool resolved = false;
    double frequency_hz = 0, harmonic_depth_db = 0, percentile_depth_db = 0;
};
// Descriptive harmonic-fit metric, not the IOA reference method or an IEC rating.
Modulation modulation(const std::vector<double> &levels, double dt, double min_hz, double max_hz);
double apparent_sound_power(double pressure_level, double distance, double reflection_correction_db = 0);
struct ArrivalSample {
    double reception, emission, wind;
    std::vector<double> power;
    double frequency_hz = 0; // used only by independent tonal sources
    double distance_m = 0;   // tone path length for absorption at received frequency
};
class ArrivalSeries {
    std::vector<ArrivalSample> samples_;

  public:
    void append(ArrivalSample);
    const std::vector<ArrivalSample> &samples() const noexcept { return samples_; }
    ArrivalSample at(double receiver_time) const;
    // dt_emission/dt_reception, valid inside the sampled subsonic trajectory.
    double doppler(double receiver_time) const;
};
class EngineeringMetrics {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    EngineeringMetrics(MetricsOptions, aeroacoustics::Parameters, std::vector<Vec3> observers,
                       std::size_t selected_nodes, Vec3 hub,
                       std::optional<aeroacoustics::PropagationOptions>);
    ~EngineeringMetrics();
    EngineeringMetrics(const EngineeringMetrics &);
    EngineeringMetrics &operator=(const EngineeringMetrics &);
    void append(double time, double hub_wind, double rotor_speed,
                const std::vector<std::vector<aeroacoustics::Node>> &, std::size_t first_node,
                std::size_t first_observer, const aeroacoustics::Snapshot &);
    void write(const std::filesystem::path &) const;
};
} // namespace turbine
