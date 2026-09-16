#include "turbine/batch.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
using namespace turbine;
namespace {
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::exception &) {
        rejected = true;
    }
    check(rejected, "Expected rejection");
}
std::string bytes(const std::filesystem::path &p) {
    std::ifstream f(p, std::ios::binary);
    check(bool(f), "Cannot read comparison file");
    return {std::istreambuf_iterator<char>(f), {}};
}
void blocks(const Case &c) {
    AcousticConfiguration config(c);
    OutputLayout layout(c, config);
    layout.blades = 1;
    layout.nodes = 5;
    layout.first = 0;
    layout.observers = 13;
    layout.parameters.freqlist = {20., 200., 2000., 16000.};
    layout.frequencies = 4;
    std::vector<aeroacoustics::Node> nodes(5);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].aero_center = {0., double(i), 80.};
        nodes[i].section.bl = {{{.001, .002}}, {{.01, .02}}, {{.002, .003}}, {{.7, 1.1}}};
        nodes[i].section.alpha_deg = -3. + i * 4.;
    }
    std::vector<Vec3> observers(13);
    for (std::size_t i = 0; i < observers.size(); ++i)
        observers[i] = {175., double(i * 10), 2.};
    for (int model : {1, 2}) {
        auto p = layout.parameters;
        p.tbltemod = model;
        aeroacoustics::AcousticWorkspace workspace(p);
        const auto full = workspace.evaluate(nodes, observers);
        layout.output_count = 4;
        const auto expected = AcousticAggregator(layout).aggregate(0, full);
        for (int outputs = 1; outputs <= 4; ++outputs) {
            layout.output_count = outputs;
            for (std::size_t size : {1, 4, 13, 20}) {
                AcousticAggregator aggregate(layout);
                aggregate.begin();
                std::size_t seen = 0;
                workspace.evaluate_blocks(
                    nodes, observers,
                    [&](std::size_t first, const auto &block) {
                        check(first == seen && block.size() <= size, "Wrong observer block");
                        for (std::size_t i = 0; i < block.size(); ++i)
                            check(block[i] == full[first + i], "Block differs from full snapshot");
                        aggregate.append(0, first, block);
                        seen += block.size();
                    },
                    size);
                const auto &actual = aggregate.finish();
                for (int k = 0; k < 4; ++k)
                    check(k < outputs ? actual.power[k] == expected.power[k] : actual.power[k].empty(),
                          "Selective aggregation changed requested output or allocated unused output");
                rejects([&] { aggregate.append(0, 0, full); });
                aggregate.begin();
                rejects([&] { aggregate.finish(); });
            }
        }
        rejects([&] { workspace.evaluate_blocks(nodes, observers, {}, 1); });
        check(workspace.evaluate(nodes, {}).empty(), "Empty observer request returned stale data");
    }
}
void solvers(Case c) {
    double worst = 0;
    std::size_t reference_calls = 0, scaled_calls = 0;
    for (double dt : {.003125, .00625, .0125})
        for (double wind : {8., 9.}) {
            c.dt = dt;
            c.wind.speed = wind;
            SolverOptions options;
            options.mode = SolverMode::scaled;
            TurbineModel model(c);
            Solver reference(model), scaled(model, options);
            options.reuse_jacobian = false;
            Solver fresh(model, options);
            const auto compare = [&](const Solver &a, const Solver &b) {
                for (int blade = 0; blade < 3; ++blade)
                    for (int j = 0; j < 3; ++j) {
                        const double dq = std::abs(a.state()[blade].q[j] - b.state()[blade].q[j]);
                        const double dv = std::abs(a.state()[blade].qd[j] - b.state()[blade].qd[j]);
                        worst = std::max({worst, dq, dv});
                        check(dq < 1e-8 && dv < 1e-8, "Scaled structural solution differs");
                    }
            };
            for (int i = 0; i < int(1. / dt); ++i) {
                reference.step();
                scaled.step();
                fresh.step();
                compare(reference, scaled);
                compare(fresh, scaled);
                check(scaled.diagnostics().scaled_residual <= 1., "Scaled residual did not converge");
            }
            reference_calls += reference.diagnostics().acceleration_evaluations;
            scaled_calls += scaled.diagnostics().acceleration_evaluations;
            auto checkpoint = scaled.checkpoint();
            scaled.step();
            auto state = scaled.state();
            scaled.restore(checkpoint);
            scaled.step();
            for (int b = 0; b < 3; ++b)
                check(state[b].q == scaled.state()[b].q, "Scaled checkpoint replay changed state");
            scaled.reset();
            check(scaled.options().mode == SolverMode::scaled, "Reset lost solver options");
            options.max_iterations = 1;
            Solver limited(model, options);
            try {
                limited.step();
                check(false, "Iteration limit not enforced");
            } catch (const std::runtime_error &e) {
                check(limited.failed() && std::string(e.what()).find("scaled_residual=") != std::string::npos,
                      "Failure lacks convergence context");
            }
        }
    check(scaled_calls < reference_calls, "Jacobian reuse did not reduce acceleration evaluations");
    SolverOptions bad;
    bad.perturbation = 0;
    rejects([&] { Solver invalid(c, bad); });
    bad = {};
    bad.acceleration_scale[0] = std::numeric_limits<double>::infinity();
    rejects([&] { Solver invalid(c, bad); });
    std::cout << "solver max_state_difference=" << worst
              << " acceleration_evaluations reference=" << reference_calls << " scaled=" << scaled_calls
              << '\n';
}
void turbulence_blocks() {
    using namespace aeroacoustics;
    for (int method : {1, 2}) {
        Parameters p;
        p.freqlist = {100., 1000.};
        Spectrum span{1., 2., 3., 4., 5.};
        std::vector<Vec3> observers{{175., 0., 2.}, {0., 175., 2.}, {-175., 0., 2.}};
        AcousticDriver full(p, span, 3, observers, .1, 0., 70., 80., method);
        AcousticDriver blocked = full;
        std::vector<std::vector<Node>> blades(3, std::vector<Node>(5));
        // More than five seconds exercises the filled moving TI window as well.
        for (int step = 0; step <= 960; ++step) {
            for (auto &blade : blades)
                for (std::size_t i = 0; i < blade.size(); ++i) {
                    blade[i].aero_center = {0., 0., 81. + i};
                    blade[i].inflow = {8. + .3 * std::sin(step * .1), 0., 0.};
                    blade[i].section.speed = 40. + i + .1 * std::sin(step * .02);
                }
            const double time = step * .00625;
            const auto *expected = full.step_view(time, blades);
            std::size_t seen = 0;
            const bool sampled = blocked.step_blocks(
                time, blades,
                [&](std::size_t first, const auto &block) {
                    check(expected && first == seen, "TI driver block order differs");
                    for (std::size_t i = 0; i < block.size(); ++i)
                        check(block[i] == (*expected)[first + i], "TI sampling order changed");
                    seen += block.size();
                },
                2);
            check(sampled == bool(expected) && seen == (sampled ? observers.size() : 0), "Sampling differs");
            check(full.turbulence_state().values == blocked.turbulence_state().values,
                  "Non-sampling TI state changed");
        }
    }
}
void batch(const std::filesystem::path &input, const std::filesystem::path &out) {
    std::vector<CaseJob> jobs;
    RunOptions options;
    options.duration = .1;
    for (int i = 0; i < 3; ++i)
        jobs.push_back({input, out / ("serial" + std::to_string(i)), options});
    const auto serial = run_cases(jobs, 1);
    for (int i = 0; i < 3; ++i)
        jobs[i].output = out / ("parallel" + std::to_string(i));
    const auto parallel = run_cases(jobs, 3);
    for (int i = 0; i < 3; ++i) {
        check(serial[i].summary.has_value() && parallel[i].summary.has_value(), "Batch case failed");
        check(serial[i].summary->lookup.calls == parallel[i].summary->lookup.calls, "Lookup sessions mixed");
        for (const auto &entry : std::filesystem::directory_iterator(out / ("serial" + std::to_string(i)))) {
            if (entry.path().filename() == "run.json")
                continue;
            check(bytes(entry.path()) == bytes(jobs[i].output / entry.path().filename()),
                  "Parallel output differs");
        }
    }
    jobs[0].input = out / "missing.fst";
    for (int i = 0; i < 3; ++i)
        jobs[i].output = out / ("mixed" + std::to_string(i));
    const auto mixed = run_cases(jobs, 2);
    check(!mixed[0].summary && !mixed[0].error.empty() && mixed[1].summary && mixed[2].summary,
          "Batch did not isolate failure");
    jobs[1].output = jobs[0].output / "nested";
    rejects([&] { run_cases(jobs, 2); });
    rejects([&] { run_cases({}, 0); });
    check(run_cases({}, 2).empty(), "Empty batch failed");
}
} // namespace
int main(int argc, char **argv) {
    try {
        check(argc == 3, "Usage: performance_api_probe CASE.fst OUTPUT");
        Case c(argv[1]);
        blocks(c);
        turbulence_blocks();
        solvers(c);
        batch(argv[1], argv[2]);
        std::cout << "P4-P7 API checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
