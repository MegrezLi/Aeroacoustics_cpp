#include "turbine/simulation.hpp"
#include "turbine/file_output.hpp"
#include <chrono>
namespace turbine {
namespace {
AcousticConfiguration configuration(const Case &c, const RunOptions &options) {
    AcousticConfiguration result(c);
    if (!options.observers.empty()) {
        if (options.metrics && !options.metrics->observers.empty())
            throw std::invalid_argument("Duplicate observer overrides");
        for (const auto &p : options.observers)
            for (double x : p)
                if (!std::isfinite(x))
                    throw std::invalid_argument("Invalid observer");
        result.observers = options.observers;
    }
    if (options.metrics && !options.metrics->observers.empty())
        result.observers = options.metrics->observers;
    return result;
}
} // namespace
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
    std::optional<EngineeringMetrics> metrics;
    double duration;
    bool started = false, finished = false, failed = false;
    std::size_t samples = 0;
    Impl(TurbineModel m, RunOptions o)
        : model(o.surfaces ? TurbineModel(o.surfaces->apply(m.data())) : std::move(m)), options(o),
          config(configuration(model.data(), o)), layout(model.data(), config),
          adapter(model.data(), config, o.surfaces), acoustic(config.make_driver()), aggregator(layout),
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
        layout.time_dependent_wind = bool(options.solver.wind);
        layout.controller = options.solver.controller;
        layout.propagation = options.propagation;
        layout.surfaces = options.surfaces;
        layout.tower = options.solver.tower;
        if (options.metrics)
            metrics.emplace(*options.metrics, config.parameters, config.observers,
                            config.blades * (config.span.size() - acoustic.first_node()),
                            solver->rotor().structure().hub, options.propagation);
        if (options.propagation)
            acoustic.set_propagation(std::make_shared<const aeroacoustics::OutdoorPropagation>(
                *options.propagation, layout.acoustic_metadata->bands));
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
                if (s.metrics) {
                    const auto &rotor = solver.rotor();
                    const auto wind = rotor.wind_at(solver.time(), rotor.structure().hub);
                    s.metrics->append(solver.time(), std::hypot(wind[0], wind[1]), rotor.structure().omega,
                                      nodes, s.acoustic.first_node(), first, block);
                }
            },
            s.options.observer_block_size);
        const AcousticResult *result = nullptr;
        if (sampled) {
            result = &s.aggregator.finish();
            ++s.samples;
        }
        s.started = true;
        return StepView{solver.time(), solver.state(), result, &solver.generalized_state(),
                        solver.operating_state()};
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
            s.solver->diagnostics(),
            s.metrics ? std::make_shared<const EngineeringMetrics>(*s.metrics) : nullptr};
}
void Simulation::reset() {
    // Reset from the loaded configuration and tables; do not reopen input files.
    auto fresh = std::make_unique<Impl>(*impl_);
    fresh->began = std::chrono::steady_clock::now();
    fresh->lookup = {};
    fresh->lookup.policy = fresh->options.lookup_policy;
    diagnostics::LookupSession session(fresh->lookup);
    fresh->acoustic = fresh->config.make_driver();
    if (fresh->options.propagation)
        fresh->acoustic.set_propagation(std::make_shared<const aeroacoustics::OutdoorPropagation>(
            *fresh->options.propagation, fresh->layout.acoustic_metadata->bands));
    fresh->solver.emplace(fresh->model, fresh->options.solver);
    if (fresh->options.metrics)
        fresh->metrics.emplace(*fresh->options.metrics, fresh->config.parameters, fresh->config.observers,
                               fresh->config.blades *
                                   (fresh->config.span.size() - fresh->acoustic.first_node()),
                               fresh->solver->rotor().structure().hub, fresh->options.propagation);
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
