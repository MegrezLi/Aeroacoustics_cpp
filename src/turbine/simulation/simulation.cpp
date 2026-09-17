#include "turbine/simulation.hpp"
#include "turbine/file_output.hpp"
#include <chrono>
namespace turbine {
struct Simulation::Impl {
    std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
    TurbineModel model;
    RunOptions options;
    diagnostics::LookupReport lookup;
    AcousticConfiguration config;
    OutputLayout layout;
    AcousticInputAdapter adapter;
    aeroacoustics::AcousticDriver acoustic;
    AcousticAggregator aggregator;
    std::optional<Solver> solver;
    double duration;
    bool started = false, finished = false, failed = false;
    std::size_t samples = 0;
    Impl(TurbineModel m, RunOptions o)
        : model(std::move(m)), options(o), config(model.data()), layout(model.data(), config),
          adapter(model.data(), config), acoustic(config.make_driver()), aggregator(layout),
          duration(o.duration.value_or(model.data().duration)) {
        if (!options.observer_block_size)
            throw std::invalid_argument("Observer block size must be positive");
        if (!std::isfinite(duration) || duration < 0)
            throw std::invalid_argument("Duration must be finite and nonnegative");
        lookup.policy = options.lookup_policy;
        diagnostics::LookupSession session(lookup);
        solver.emplace(model, options.solver);
        layout.dofs = solver->layout().dofs();
        layout.coupling_blocks = solver->layout().blocks();
    }
};
Simulation::Simulation(TurbineModel m, RunOptions o) : impl_(std::make_unique<Impl>(std::move(m), o)) {}
Simulation::~Simulation() = default;
Simulation::Simulation(const Simulation &other) : impl_(std::make_unique<Impl>(*other.impl_)) {}
Simulation &Simulation::operator=(const Simulation &other) {
    if (this != &other)
        impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}
Simulation::Simulation(Simulation &&) noexcept = default;
Simulation &Simulation::operator=(Simulation &&) noexcept = default;
const OutputLayout &Simulation::layout() const { return impl_->layout; }
const Solver &Simulation::solver() const { return *impl_->solver; }
std::optional<StepView> Simulation::next() {
    auto &s = *impl_;
    if (s.failed)
        throw std::logic_error("Simulation failed; reset or restore before advancing");
    if (s.finished)
        return std::nullopt;
    diagnostics::LookupSession session(s.lookup);
    try {
        auto &solver = *s.solver;
        if (s.started) {
            if (solver.time() + s.model.data().dt / 2 >= s.duration) {
                s.finished = true;
                return std::nullopt;
            }
            solver.step();
        }
        const auto &nodes =
            s.adapter.update(solver.aerodynamic(), solver.time(), s.acoustic.is_sample_time(solver.time()),
                             s.acoustic.first_node());
        if (s.acoustic.is_sample_time(solver.time()))
            s.aggregator.begin();
        const bool sampled = s.acoustic.step_blocks(
            solver.time(), nodes,
            [&](std::size_t first, const aeroacoustics::Snapshot &block) {
                s.aggregator.append(solver.time(), first, block);
            },
            s.options.observer_block_size);
        const AcousticResult *result = nullptr;
        if (sampled) {
            result = &s.aggregator.finish();
            ++s.samples;
        }
        s.started = true;
        return StepView{solver.time(), solver.state(), result, &solver.generalized_state()};
    } catch (...) {
        s.failed = true;
        throw;
    }
}
RunSummary Simulation::summary() const {
    const auto &s = *impl_;
    if (!s.finished || s.failed)
        throw std::logic_error("Simulation has not completed successfully");
    return {s.model.data().dt,
            s.solver->time(),
            std::chrono::duration<double>(std::chrono::steady_clock::now() - s.began).count(),
            s.solver->step_number(),
            s.samples,
            s.lookup,
            s.solver->diagnostics()};
}
void Simulation::reset() {
    // Reset from the loaded configuration and tables; do not reopen input files.
    auto fresh = std::make_unique<Impl>(*impl_);
    fresh->began = std::chrono::steady_clock::now();
    fresh->lookup = {};
    fresh->lookup.policy = fresh->options.lookup_policy;
    diagnostics::LookupSession session(fresh->lookup);
    fresh->acoustic = fresh->config.make_driver();
    fresh->solver.emplace(fresh->model, fresh->options.solver);
    fresh->started = fresh->finished = fresh->failed = false;
    fresh->samples = 0;
    impl_ = std::move(fresh);
}
Simulation::Checkpoint Simulation::checkpoint() const {
    if (impl_->failed)
        throw std::logic_error("Cannot checkpoint a failed simulation");
    return Checkpoint(std::make_shared<const Impl>(*impl_));
}
void Simulation::restore(const Checkpoint &snapshot) {
    if (!impl_->model.same_model(snapshot.saved_->model))
        throw std::invalid_argument("Checkpoint belongs to a different turbine model");
    impl_ = std::make_unique<Impl>(*snapshot.saved_);
}
RunSummary run(Simulation &simulation, ResultSink &sink) {
    if (simulation.impl_->started || simulation.impl_->failed)
        throw std::logic_error("run() requires a fresh/reset simulation; use next() for resumed segments");
    sink.begin(simulation.layout());
    while (auto frame = simulation.next())
        sink.write(*frame);
    auto result = simulation.summary();
    sink.finish(result);
    return result;
}
RunSummary run_case(const std::filesystem::path &input, const std::filesystem::path &output,
                    RunOptions options) {
    TurbineModel model{Case(input)};
    FileOutput sink(output);
    Simulation simulation(std::move(model), options);
    return run(simulation, sink);
}
} // namespace turbine
