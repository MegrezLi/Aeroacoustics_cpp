#pragma once
#include "checked_output.hpp"
#include "turbine/results.hpp"
namespace turbine {
class FileOutput final : public ResultSink {
    std::filesystem::path directory_;
    diagnostics::CheckedOutput metadata_, dynamics_, lookup_output_;
    std::array<diagnostics::CheckedOutput, 4> outputs_, masks_;
    std::optional<OutputLayout> layout_;
    bool finished_ = false, failed_ = false;
    double last_time_ = -std::numeric_limits<double>::infinity();
    void begin_output(const OutputLayout &);
    void write_step(const StepView &);
    void finish_output(const RunSummary &);

  public:
    // Truncates run.json immediately; success is recorded only after all writes finish.
    explicit FileOutput(std::filesystem::path directory);
    void begin(const OutputLayout &) override;
    void write(const StepView &) override;
    void finish(const RunSummary &) override;
};
} // namespace turbine
