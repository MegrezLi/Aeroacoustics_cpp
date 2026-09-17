#include "turbine/simulation.hpp"
#include <iostream>
#include <type_traits>
using namespace aeroacoustics;
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
void near(double a, double b, double tolerance = 1e-12) {
    check(std::abs(a - b) <= tolerance, "Numerical semantic check failed");
}
void quantities() {
    static_assert(!std::is_convertible_v<PressurePsd, BandMeanSquarePressure>);
    static_assert(!std::is_convertible_v<BandSoundPower, BandMeanSquarePressure>);
    static_assert(!std::is_convertible_v<BandSoundPowerLevel, BandSoundPressureLevel>);
    const FrequencyBands bands({{100, 80, 120}, {150, 120, 180}});
    const auto pressure = integrate(PressurePsd(bands, {2e-6, 2e-6}));
    near(pressure.values()[0], 8e-5);
    near(total_mean_square_pressure(pressure), 2e-4);
    auto level = pressure_levels(BandMeanSquarePressure(bands, {4e-10, 0}));
    near(level.values()[0], 0);
    check(level.values()[1] == -INFINITY, "Silence must remain negative infinity");
    const auto restored = mean_square_pressure(level);
    near(restored.values()[0], 4e-10, 1e-24);
    near(restored.values()[1], 0);
    near(power_levels(BandSoundPower(bands, {1e-12, 1.})).values()[1], 120.);
    near(sound_power(power_levels(BandSoundPower(bands, {1e-12, 1.}))).values()[1], 1.);
    check(std::isfinite(mean_square_pressure(BandSoundPressureLevel(bands, {3100, 0})).values()[0]),
          "Finite pressure overflowed during reference conversion");
    const auto third = FrequencyBands::third_octave(-1, 3);
    const FrequencyBands octave({{1000, third.values().front().lower_hz, third.values().back().upper_hz}});
    const auto merged = merge_bands(BandMeanSquarePressure(third, {1, 2, 3}), octave);
    near(merged.values()[0], 6.);
    near(total_sound_power(merge_bands(BandSoundPower(third, {1, 2, 3}), octave)), 6.);
    const auto a = apply_a_weighting(BandMeanSquarePressure(third, {1, 2, 3}));
    const auto weights = a_weighting({third.values()[0].center_hz, 1000., third.values()[2].center_hz});
    for (std::size_t i = 0; i < 3; ++i)
        near(a.values()[i], (i + 1) * std::pow(10., weights[i] / 10.));
    rejects([&] { apply_a_weighting(a); });
    rejects([&] { FrequencyBands invalid({{100, 80, 120}, {110, 100, 130}}); });
    rejects([&] { PressurePsd invalid(bands, {1}); });
    rejects([&] { PressurePsd invalid(bands, {-1, 1}); });
    rejects([&] { merge_bands(pressure, FrequencyBands({{120, 90, 180}})); });
    rejects([&] {
        merge_bands(BandMeanSquarePressure(FrequencyBands({{100, 80, 120}, {160, 140, 180}}), {1, 2}),
                    FrequencyBands({{120, 80, 180}}));
    });
    Parameters p;
    p.freqlist = {100., 1000., 10000.};
    const auto typed = band_section_spectrum(p, Section{});
    const auto raw = section_spectrum(p, Section{});
    for (std::size_t m = 0; m < raw.size(); ++m)
        check(typed[m].values() == raw[m], "Typed wrapper changed kernel values");
    p.freqlist = {100., 101.};
    section_spectrum(p, Section{});
    rejects([&] { band_section_spectrum(p, Section{}); });
    rejects([&] { FrequencyBands::openfast_reference(p.freqlist); });
    for (double f : openfast_centers_hz) {
        const double ratio = std::pow(2., 1. / 3.), omega = 2 * 3.14159265358979323846 * f;
        check(reference_tno_bandwidth(f).value == 2 * omega * (std::sqrt(ratio) - 1 / std::sqrt(ratio)),
              "TNO reference convention changed");
    }
}
struct CoupledOscillator : AccelerationOperator {
    void evaluate(std::size_t block, StateView s, double *a) const override {
        double sum = 0;
        for (std::size_t j = 0; j < s.size; ++j)
            sum += s.q[j];
        for (std::size_t i = 0; i < s.size; ++i)
            a[i] = double(i + 1 + block) - 2 * s.q[i] - .3 * sum - .2 * s.qd[i];
    }
};
DofLayout layout(std::size_t n) {
    std::vector<DofDescriptor> dofs;
    for (std::size_t i = 0; i < n; ++i)
        dofs.push_back({"coordinate" + std::to_string(i), "m", 1.});
    return {std::move(dofs), {{"fully_coupled", 0, n}}};
}
void integrators() {
    CoupledOscillator op;
    for (std::size_t n : {1, 2, 3, 4, 6}) {
        GeneralizedAlpha solver(layout(n), .01, .5);
        rejects([&] { solver.correct(op); });
        solver.predict();
        solver.correct(op);
        // First step has a zero predictor. Invert the diagonal + rank-one
        // Jacobian analytically (independent of the implementation's dense solve).
        const double am = 0., af = 1. / 3., gamma = .5 - am + af, beta = std::pow(1 - am + af, 2) / 4;
        const double bp = .01 * .01 * beta * (1 - af) / (1 - am), gp = .01 * gamma * (1 - af) / (1 - am);
        const double d = 1 + 2 * bp + .2 * gp, c = .3 * bp, sum = n * (n + 1) / 2.;
        for (std::size_t i = 0; i < n; ++i) {
            const double a = (i + 1) / d - c * sum / (d * (d + c * n));
            near(solver.acceleration()[i], a, 1e-11);
            near(solver.state().q[i], bp * a);
            near(solver.state().qd[i], gp * a);
        }
        for (int step = 0; step < 15; ++step) {
            solver.predict();
            solver.correct(op);
        }
        auto copy = solver;
        solver.predict();
        solver.correct(op);
        copy.predict();
        copy.correct(op);
        check(solver.state().q == copy.state().q && solver.state().qd == copy.state().qd,
              "Integrator copy lost history");
        copy.reset();
        near(copy.time(), 0);
        check(copy.step_number() == 0, "Reset retained step count");
        for (double x : copy.state().q)
            near(x, 0);
        copy.predict();
        rejects([&] { copy.predict(); });
    }
    rejects([&] { DofLayout invalid({{"x", "m", 1}, {"y", "m", 1}}, {{"a", 0, 1}}); });
    rejects([&] { DofLayout invalid({{"x", "m", 1}, {"x", "m", 1}}, {{"a", 0, 2}}); });
    rejects([&] { DofLayout invalid({{"x", "m", 0}}, {{"a", 0, 1}}); });
    rejects([&] { GeneralizedAlpha invalid(layout(4), 0, .5); });
    struct PartialOutput : AccelerationOperator {
        void evaluate(std::size_t, StateView, double *a) const override { a[0] = 0; }
    } partial;
    GeneralizedAlpha failed(layout(4), .01, .5);
    failed.predict();
    rejects([&] { failed.correct(partial); });
    check(failed.failed(), "Partial/non-finite module output not rejected");
    rejects([&] { failed.predict(); });
    failed.reset();
    failed.predict();
    failed.correct(op);
}
void configuration(const std::filesystem::path &path) {
    Case c(path);
    check(configure_modules(c).blades == 3, "Wrong module profile");
    check(module_capabilities().size() == 6, "Missing module contracts");
    auto incompatible = configure_modules(c);
    incompatible.modules[0] = incompatible.modules[1];
    rejects([&] { structural_dof_layout(incompatible); });
    auto bad = c;
    bad.structure.rows[bad.structure.line("TwFADOF1")][0] = "True";
    rejects([&] { configure_modules(bad); });
    rejects([&] { TurbineModel invalid(bad); });
    bad = c;
    bad.primary.rows[bad.primary.line("CompServo")][0] = "1";
    rejects([&] { configure_modules(bad); });
    const auto dofs = fixed_base_dof_layout();
    check(dofs.size() == 9 && dofs.blocks().size() == 3, "Wrong fixed-base DOF layout");
    Simulation simulation(c, RunOptions{.1});
    while (auto frame = simulation.next()) {
        check(frame->generalized && frame->generalized->q.size() == 9, "Missing generalized state view");
        for (std::size_t b = 0; b < 3; ++b)
            for (std::size_t j = 0; j < 3; ++j)
                check(frame->state[b].q[j] == frame->generalized->q[3 * b + j], "Legacy state view differs");
        if (frame->acoustics)
            check(frame->acoustics->metadata && frame->acoustics->metadata->bands.size() == 34,
                  "Missing acoustic metadata");
    }
    AcousticConfiguration acoustic(c);
    OutputLayout out(c, acoustic);
    out.parameters.freqlist = {100, 101};
    out.frequencies = 2;
    rejects([&] { AcousticAggregator invalid(out); });
}
} // namespace
int main(int argc, char **argv) {
    try {
        check(argc == 2, "Usage: semantics_modules_probe CASE.fst");
        quantities();
        integrators();
        configuration(argv[1]);
        std::cout << "S4-S5 quantities, layouts, module contracts and integrator checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
