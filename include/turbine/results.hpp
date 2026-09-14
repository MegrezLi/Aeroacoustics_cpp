#pragma once
#include "lookup_diagnostics.hpp"
#include "turbine/acoustic_adapter.hpp"
namespace turbine {
struct OutputLayout {
    aeroacoustics::Parameters parameters;
    std::array<std::vector<std::string>, 4> labels;
    std::string prefix;
    int output_count;
    std::size_t blades, nodes, first, observers, frequencies;
    OutputLayout(const Case &, const AcousticConfiguration &);
};
struct AcousticResult {
    // Relative mean-square pressure, not dB. Order: total, bands, mechanisms, nodes.
    std::array<std::vector<double>, 4> power;
};
class AcousticAggregator {
    OutputLayout layout_;
    AcousticResult result_;

  public:
    explicit AcousticAggregator(OutputLayout);
    const AcousticResult &aggregate(double time, const aeroacoustics::Snapshot &);
};
struct StepView {
    double time;
    const RotorState &state;
    const AcousticResult *acoustics; // null on non-sampling steps
};
struct RunSummary {
    double dt, duration, elapsed_seconds;
    std::size_t steps, acoustic_samples;
    diagnostics::LookupReport lookup;
};
class ResultSink {
  public:
    virtual ~ResultSink() = default;
    virtual void begin(const OutputLayout &) = 0;
    virtual void write(const StepView &) = 0;
    virtual void finish(const RunSummary &) = 0;
};
} // namespace turbine
