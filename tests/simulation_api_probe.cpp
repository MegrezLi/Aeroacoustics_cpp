#include "acoustic_levels.hpp"
#include "turbine/file_output.hpp"
#include "turbine/simulation.hpp"
#include <fstream>
#include <iostream>
#include <type_traits>
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
struct MemorySink : ResultSink {
    std::vector<double> values;
    std::size_t frames = 0, samples = 0;
    bool finished = false;
    void begin(const OutputLayout &layout) override {
        check(layout.labels[2].size() ==
                  layout.observers * layout.frequencies * aeroacoustics::mechanism_count,
              "Mechanism channel count");
        values.clear();
        frames = samples = 0;
        finished = false;
    }
    void write(const StepView &frame) override {
        values.push_back(frame.time);
        ++frames;
        for (const auto &blade : frame.state)
            for (const auto &v : {blade.q, blade.qd})
                values.insert(values.end(), v.begin(), v.end());
        values.push_back(frame.acoustics ? 1 : 0);
        if (frame.acoustics) {
            ++samples;
            for (const auto &v : frame.acoustics->power)
                values.insert(values.end(), v.begin(), v.end());
        }
    }
    void finish(const RunSummary &) override { finished = true; }
};
std::vector<double> signature(const Solver &s) {
    std::vector<double> out{s.time(), double(s.step_number())};
    for (const auto &b : s.state())
        for (const auto &v : {b.q, b.qd})
            out.insert(out.end(), v.begin(), v.end());
    for (const auto &blade : s.aerodynamic().blades)
        for (const auto &a : blade) {
            for (double x : {a.phi, a.alpha, a.speed, a.axial, a.tangential, a.coefficients.cl,
                             a.coefficients.cd, a.coefficients.cm})
                out.push_back(x);
            for (const auto &v : {a.motion.position, a.motion.velocity, a.load.force, a.load.moment})
                out.insert(out.end(), v.begin(), v.end());
        }
    return out;
}
void drain(Simulation &s, MemorySink &sink) {
    while (auto frame = s.next())
        sink.write(*frame);
}
} // namespace
static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<Solver &>().state())>>);
static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<Rotor &>().structure())>>);
int main(int argc, char **argv) {
    try {
        check(argc == 3, "Usage: simulation_api_probe CASE.fst OUTPUT_DIR");
        TurbineModel model{Case(argv[1])};
        Solver reference(model);
        auto isolated = [&] {
            Case original(argv[1]);
            Solver result(original);
            original.airfoils.clear();
            original.stations.clear();
            original.wind.speed = 99;
            return result;
        }();
        for (int i = 0; i < 12; ++i) {
            reference.step();
            isolated.step();
        }
        check(signature(reference) == signature(isolated), "Solver retained mutable external Case");
        auto snapshot = reference.checkpoint();
        auto copied = reference;
        for (int i = 0; i < 20; ++i) {
            reference.step();
            copied.step();
        }
        const auto expected = signature(reference);
        check(signature(copied) == expected, "Solver copy shares evolving state");
        copied.restore(snapshot);
        for (int i = 0; i < 20; ++i)
            copied.step();
        check(signature(copied) == expected, "Solver checkpoint lost UA/induction history");
        rejects([&] { isolated.restore(snapshot); });
        copied.reset();
        check(signature(copied) == signature(Solver(model)), "Solver reset differs from initial state");
        Solver moved(std::move(copied));
        moved.step();
        Solver first(model);
        first.step();
        check(signature(moved) == signature(first), "Solver move lost owned model");

        auto ua = [&] {
            auto foil = model.data().airfoils.at(10);
            UnsteadyAirfoil result(foil, 2, .00625, 340);
            foil.alpha.clear();
            foil.coefficients.clear();
            return result;
        }();
        UnsteadyAirfoil ua_reference(model.data().airfoils.at(10), 2, .00625, 340);
        for (int i = 0; i < 20; ++i) {
            const double a = .05 + .001 * i;
            ua.advance(a, 50, i);
            ua_reference.advance(a, 50, i);
            auto x = ua.evaluate(a, 50), y = ua_reference.evaluate(a, 50);
            check(x.cl == y.cl && x.cd == y.cd && x.cm == y.cm, "UA retained destroyed Airfoil");
        }

        RunOptions options;
        options.duration = .4;
        Simulation full(model, options);
        MemorySink all;
        const auto summary = run(full, all);
        check(all.frames == 65 && all.samples == 5 && all.finished, "Library time/sample schedule");
        check(summary.steps == 64 && summary.acoustic_samples == 5, "Library summary");
        Simulation segment(model, options);
        MemorySink prefix;
        for (int i = 0; i < 19; ++i)
            prefix.write(*segment.next());
        auto saved = segment.checkpoint();
        Simulation branch = segment;
        MemorySink suffix;
        drain(segment, suffix);
        prefix.values.insert(prefix.values.end(), suffix.values.begin(), suffix.values.end());
        check(prefix.values == all.values, "Library segmented output differs");
        MemorySink copied_suffix;
        drain(branch, copied_suffix);
        check(copied_suffix.values == suffix.values, "Simulation copy lost TI/UA state");
        segment.restore(saved);
        MemorySink replay;
        drain(segment, replay);
        check(replay.values == suffix.values, "Simulation checkpoint replay differs");
        segment.reset();
        MemorySink reset;
        run(segment, reset);
        check(reset.values == all.values, "Simulation reset differs");
        Simulation foreign(Case(argv[1]), options);
        rejects([&] { foreign.restore(saved); });
        Simulation moved_sim(std::move(branch));
        check(!moved_sim.next(), "Moved completed simulation changed");
        Simulation partial(model, options);
        partial.next();
        rejects([&] {
            MemorySink sink;
            run(partial, sink);
        });
        Case supersonic_input(argv[1]);
        supersonic_input.sound_speed = 1;
        Solver failed_solver(supersonic_input);
        failed_solver.step(); // UA deliberately skips its history advance at step zero.
        auto solver_before_failure = failed_solver.checkpoint();
        rejects([&] { failed_solver.step(); });
        check(failed_solver.failed(), "Solver did not record failure");
        rejects([&] { failed_solver.checkpoint(); });
        failed_solver.restore(solver_before_failure);
        check(!failed_solver.failed(), "Solver restore did not clear failure");
        Simulation failed_simulation(supersonic_input, options);
        auto before_failure = failed_simulation.checkpoint();
        rejects([&] { failed_simulation.next(); }); // Supersonic acoustic input.
        rejects([&] { failed_simulation.next(); });
        rejects([&] { failed_simulation.checkpoint(); });
        failed_simulation.restore(before_failure);
        failed_simulation.checkpoint();
        failed_simulation.reset();
        failed_simulation.checkpoint();

        // Independent diagnostic sessions remain intact when simulations interleave.
        diagnostics::LookupReport outer;
        diagnostics::LookupSession outer_session(outer);
        Simulation one(model, options), two(model, options);
        for (int i = 0; i < 4; ++i) {
            one.next();
            two.next();
        }
        check(diagnostics::active_report == &outer && outer.calls == 0, "Simulation leaked lookup context");

        const std::filesystem::path output = argv[2];
        FileOutput files(output);
        Simulation written(model, options);
        run(written, files);
        check(std::filesystem::file_size(output / "run.json") > 0, "File sink did not finish");
        FileOutput bad(output / "bad");
        bad.begin(written.layout());
        AcousticResult invalid;
        invalid.power[0] = {0};
        StepView invalid_frame{0, written.solver().state(), &invalid};
        rejects([&] { bad.write(invalid_frame); });
        rejects([&] { bad.finish(summary); });
        check(std::filesystem::file_size(output / "bad" / "run.json") == 0,
              "Failed sink wrote success metadata");
        const auto input_copy = output / "owned-inputs";
        const auto moved_inputs = output / "moved-inputs";
        check(!std::filesystem::exists(input_copy) && !std::filesystem::exists(moved_inputs),
              "Ownership test requires fresh input-copy directories");
        std::filesystem::copy(std::filesystem::path(argv[1]).parent_path(), input_copy,
                              std::filesystem::copy_options::recursive);
        Simulation detached(Case(input_copy / std::filesystem::path(argv[1]).filename()), options);
        detached.next();
        // Both names are fixed children of the test output; official inputs remain untouched.
        std::filesystem::rename(input_copy, moved_inputs);
        detached.reset();
        MemorySink detached_output;
        run(detached, detached_output);
        check(detached_output.values == all.values, "Reset reopened input files instead of owned data");
        std::cout << "S1-S3 simulation, ownership, checkpoint and output checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
