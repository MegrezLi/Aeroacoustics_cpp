#pragma once
#include "turbine/simulation.hpp"
namespace turbine {
struct CaseJob {
    std::filesystem::path input, output;
    RunOptions options;
};
struct CaseOutcome {
    std::optional<RunSummary> summary;
    std::string error;
};
// Results follow input order. One failed case does not cancel other cases.
// Output directories must be distinct and non-nested; validated before any output is opened.
// Each worker owns its complete simulation, lookup session and acoustic workspace.
std::vector<CaseOutcome> run_cases(const std::vector<CaseJob> &, std::size_t workers = 1);
} // namespace turbine
