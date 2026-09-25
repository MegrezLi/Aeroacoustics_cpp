#pragma once
#include "turbine/simulation.hpp"
namespace turbine {
struct WakeSample {
    double time, inflow_speed, ct;
};
struct WakeHistory {
    Vec3 hub{};
    double radius = 1;
    std::vector<WakeSample> samples;
    void validate() const;
    WakeSample at(double time) const;
};
// Quasi-steady Jensen top-hat, aligned fixed wind; dimensional deficits use source inflow.
double jensen_deficit(double ct, double source_speed, double radius, double expansion, double downstream,
                      double transverse);
class FarmWind final : public WindField {
    SteadyWind ambient_;
    Vec3 origin_, direction_;
    double expansion_;
    std::vector<std::shared_ptr<const WakeHistory>> upstream_;

  public:
    FarmWind(SteadyWind, Vec3 origin, double expansion, std::vector<std::shared_ptr<const WakeHistory>>);
    Vec3 at(double time, const Vec3 &local_position) const override;
    const char *description() const noexcept override {
        return "fixed-direction quasi-steady Jensen/RSS farm wind";
    }
};
struct StationarySource {
    std::string name, provenance, kind; // bands or tones; Lw referenced to 1 pW
    Vec3 position{}, axis{1, 0, 0};
    double directivity = 0; // Q(theta)=1+directivity*cos(theta), integral 4*pi
    std::vector<double> frequencies, levels_db;
    void validate() const;
    static StationarySource read(const std::filesystem::path &);
    // A-weighted relative mean-square pressure, exact tone frequency for A and absorption.
    std::vector<double> receiver_power(const Vec3 &, const aeroacoustics::FrequencyBands &, double rho,
                                       double sound_speed,
                                       const std::optional<aeroacoustics::PropagationOptions> &) const;
};
struct FarmUnit {
    std::string name;
    std::filesystem::path case_file;
    Vec3 origin{}; // flat site, z must be zero; all machines retain their own hub heights
    std::optional<ControlConfig> controller;
    std::shared_ptr<const TowerInfluence> tower;
    std::shared_ptr<const SurfaceSet> surfaces;
};
struct FarmOptions {
    std::string provenance;
    std::vector<FarmUnit> turbines;
    std::vector<Vec3> observers;
    std::vector<StationarySource> sources;
    double duration = 20, statistics_start = 2, expansion = .05, max_ct = .95;
    bool wakes = true;
    std::size_t max_values = 32000000;
    std::optional<aeroacoustics::PropagationOptions> propagation;
    static FarmOptions read(const std::filesystem::path &);
};
void run_farm(const FarmOptions &, const std::filesystem::path &output);
} // namespace turbine
