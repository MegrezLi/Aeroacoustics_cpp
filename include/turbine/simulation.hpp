#pragma once
#include "propagation.hpp"
#include "turbine/results.hpp"
#include "turbine/solver.hpp"
namespace turbine {
struct RunOptions {
    std::optional<double> duration;
    diagnostics::LookupPolicy lookup_policy = diagnostics::LookupPolicy::clamp;
    std::size_t observer_block_size = 1;
    SolverOptions solver;
    std::optional<aeroacoustics::PropagationOptions> propagation;
};
class Simulation {
    friend RunSummary run(Simulation &, ResultSink &);
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    explicit Simulation(TurbineModel, RunOptions = {});
    explicit Simulation(const Case &c, RunOptions options = {}) : Simulation(TurbineModel(c), options) {}
    ~Simulation();
    Simulation(const Simulation &);
    Simulation &operator=(const Simulation &);
    Simulation(Simulation &&) noexcept;
    Simulation &operator=(Simulation &&) noexcept;
    const OutputLayout &layout() const;
    const Solver &solver() const;
    // All references expire at next(), reset(), restore(), assignment or destruction.
    std::optional<StepView> next();
    RunSummary summary() const;
    void reset();
    class Checkpoint {
        friend class Simulation;
        std::shared_ptr<const Impl> saved_;
        explicit Checkpoint(std::shared_ptr<const Impl> p) : saved_(std::move(p)) {}

      public:
        Checkpoint(const Checkpoint &) = default;
        Checkpoint &operator=(const Checkpoint &) = default;
    };
    Checkpoint checkpoint() const;
    void restore(const Checkpoint &);
};
RunSummary run(Simulation &, ResultSink &);
RunSummary run_case(const std::filesystem::path &input, const std::filesystem::path &output, RunOptions = {});
} // namespace turbine
