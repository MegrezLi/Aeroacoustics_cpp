#include "turbine/simulation.hpp"
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
using namespace turbine;
namespace {
void check(bool b, const char *m) {
    if (!b)
        throw std::runtime_error(m);
}
void near(double a, double b, double tol = 1e-10) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tol)
        throw std::runtime_error("Mismatch: " + std::to_string(a) + " vs " + std::to_string(b));
}
template <class F> void rejects(F f) {
    bool yes = false;
    try {
        f();
    } catch (const std::exception &) {
        yes = true;
    }
    check(yes, "Expected rejection");
}
std::string contents(const std::filesystem::path &p) {
    std::ifstream f(p);
    return {std::istreambuf_iterator<char>(f), {}};
}
void analytic() {
    auto s = level_statistics({0, 1, 3}, {1, 3, 3});
    near(s.duration, 3);
    near(s.leq_db, 10 * std::log10(8. / 3));
    near(s.l5_db, 10 * std::log10(3));
    near(s.l50_db, 10 * std::log10(3));
    near(s.l95_db, 10 * std::log10(1.3));
    check(level_statistics({0, 1}, {0, 0}).leq_db == -INFINITY, "Silence not preserved");
    rejects([] { level_statistics({0, 0}, {1, 2}); });
    rejects([] { level_statistics({0, 1}, {-1, 2}); });
    std::vector<double> wave;
    for (int i = 0; i < 100; ++i)
        wave.push_back(40 + 3 * std::cos(2 * pi * .5 * i * .1));
    auto am = modulation(wave, .1, .2, 1);
    check(am.resolved, "AM not resolved");
    near(am.frequency_hz, .5);
    near(am.harmonic_depth_db, 6);
    near(modulation(std::vector<double>(100, 40), .1, .2, 1).harmonic_depth_db, 0);
    rejects([&] { modulation(wave, 1, .2, 1); });
    near(apparent_sound_power(50, 10), 50 + 10 * std::log10(400 * pi));
    near(apparent_sound_power(56, 10, 6), apparent_sound_power(50, 10));
    ArrivalSeries arrival;
    // Constant approaching source: t_r=t_e+(100-10*t_e)/340.
    for (int i = 0; i < 11; ++i)
        arrival.append({i + (100. - 10 * i) / 340., double(i), 8., {2. + i}, 1000.});
    double tr = 5 + 50. / 340;
    near(arrival.at(tr).power[0], 7);
    near(arrival.doppler(tr), 340. / 330);
    rejects([&] { arrival.at(0); });
    rejects([&] { arrival.append({1, 11, 8, {1}}); });
    ArrivalSeries stationary;
    stationary.append({1, 0, 8, {4}});
    stationary.append({2, 1, 8, {4}});
    near(stationary.doppler(1.5), 1);
    near(stationary.at(1.5).power[0], 4);
}
void integration(const std::filesystem::path &out) {
    MetricsOptions o;
    o.am_max_hz = .4; // the analytic fixture below has one-second source samples
    o.start = 0;
    o.end = 1;
    o.retarded_time = false;
    aeroacoustics::Parameters p;
    p.freqlist = {1000};
    p.aweighting = true;
    EngineeringMetrics metrics(o, p, {{10, 0, 0}}, 1, {}, {});
    std::vector<std::vector<aeroacoustics::Node>> nodes(1, std::vector<aeroacoustics::Node>(1));
    aeroacoustics::Snapshot spectrum(1, std::vector<aeroacoustics::Mechanisms>(1));
    for (auto &m : spectrum[0][0])
        m = {-INFINITY};
    spectrum[0][0][0][0] = 0;
    metrics.append(0, 7.5, 1, nodes, 0, 0, spectrum);
    spectrum[0][0][0][0] = 10 * std::log10(3);
    metrics.append(1, 8.5, 1, nodes, 0, 0, spectrum);
    std::filesystem::create_directories(out);
    metrics.write(out);
    std::ifstream f(out / "wind_bins.csv");
    std::string line;
    std::getline(f, line);
    for (int bin = 0; bin < 2; ++bin) {
        check(bool(std::getline(f, line)), "Missing crossing wind bin");
        std::istringstream row(line);
        std::vector<double> values;
        std::string word;
        while (std::getline(row, word, ','))
            values.push_back(std::stod(word));
        near(values[1], 7 + bin);
        near(values[3], .5);
        near(values[4], 10 * std::log10(1.5 + bin));
    }
}
void moving_absorption(const std::filesystem::path &out) {
    MetricsOptions o;
    o.start = .2;
    o.end = .8;
    o.am_max_hz = .4;
    o.tones.push_back({"analytic", "prescribed test", 0, 0, 1000, 0, 80, 1});
    aeroacoustics::Parameters p;
    p.freqlist = {1000};
    aeroacoustics::PropagationOptions prop;
    prop.sound_speed = p.spdsound;
    EngineeringMetrics free(o, p, {{10, 0, 0}}, 1, {}, {}), absorbed(o, p, {{10, 0, 0}}, 1, {}, prop);
    std::vector<std::vector<aeroacoustics::Node>> nodes(1, std::vector<aeroacoustics::Node>(1));
    aeroacoustics::Snapshot spectrum(1, std::vector<aeroacoustics::Mechanisms>(1));
    for (auto &m : spectrum[0][0])
        m = {-INFINITY};
    for (int t = 0; t <= 1; ++t) {
        nodes[0][0].aero_center[0] = t;
        free.append(t, 8, 1, nodes, 0, 0, spectrum);
        absorbed.append(t, 8, 1, nodes, 0, 0, spectrum);
    }
    std::filesystem::create_directories(out / "free");
    std::filesystem::create_directories(out / "absorbed");
    free.write(out / "free");
    absorbed.write(out / "absorbed");
    std::ifstream a(out / "free/receiver_tones.csv"), b(out / "absorbed/receiver_tones.csv");
    std::string x, y;
    std::getline(a, x);
    std::getline(b, y);
    while (std::getline(a, x)) {
        check(bool(std::getline(b, y)), "Missing absorbed tone");
        auto parse = [](const std::string &s) {
            std::istringstream row(s);
            std::string word;
            std::vector<double> v;
            while (std::getline(row, word, ','))
                v.push_back(std::stod(word));
            return v;
        };
        auto plain = parse(x), loss = parse(y);
        double t = plain[2], factor = p.spdsound / (p.spdsound - 1), te = (t - 10 / p.spdsound) * factor;
        near(plain[3], 1000 * factor);
        near(plain[4] - loss[4], prop.atmosphere.absorption_db_per_m(plain[3]) * (10 - te));
    }
    // A failed history cannot be exported as a complete result.
    rejects([&] { free.append(0, 8, 1, nodes, 0, 0, spectrum); });
    rejects([&] { free.write(out / "free"); });
}
void surfaces(const std::filesystem::path &file, const Case &c) {
    auto set = SurfaceSet::read(file);
    const auto &d = set->data().at(0);
    auto b = d.at(0, (1000. + 100000000.) / 2, 2.);
    near(b.dstar[0], .011);
    near(b.dstar[1], .007);
    near(b.edge_velocity_ratio[0], 1.15);
    rejects([&] { d.at(171, 1e6, 1); });
    rejects([&] { d.at(0, 999, 1); });
    auto modified = set->apply(c);
    near(modified.airfoils[29].coefficients[1].cl, c.airfoils[29].coefficients[1].cl);
    aeroacoustics::Parameters p;
    p.timod = 0;
    p.bluntmod = 0;
    aeroacoustics::Section input;
    input.bl = b;
    input.tabulated_boundary_layer = true;
    auto overridden = aeroacoustics::section_spectrum(p, input);
    p.x_blmethod = 2;
    input.tabulated_boundary_layer = false;
    auto tabulated = aeroacoustics::section_spectrum(p, input);
    check(overridden == tabulated, "Per-section BL override differs from table model");
}
void simulation(const Case &c, const std::filesystem::path &surface, const std::filesystem::path &out) {
    RunOptions o;
    o.duration = 1.;
    o.surfaces = SurfaceSet::read(surface);
    o.metrics = MetricsOptions{};
    o.metrics->start = 0;
    o.metrics->end = 1;
    o.metrics->retarded_time = false;
    Simulation sim(c, o);
    sim.next();
    sim.next();
    auto cp = sim.checkpoint();
    while (sim.next()) {
    };
    auto report = sim.summary();
    std::filesystem::create_directories(out / "first");
    report.metrics->write(out / "first");
    sim.restore(cp);
    while (sim.next()) {
    };
    std::filesystem::create_directories(out / "restored");
    sim.summary().metrics->write(out / "restored");
    sim.reset();
    while (sim.next()) {
    };
    std::filesystem::create_directories(out / "reset");
    sim.summary().metrics->write(out / "reset");
    check(contents(out / "first/receiver_history.csv") == contents(out / "restored/receiver_history.csv"),
          "Metric checkpoint mismatch");
    check(contents(out / "first/receiver_history.csv") == contents(out / "reset/receiver_history.csv"),
          "Metric reset mismatch");
    auto worker = [&] {
        Simulation s(c, o);
        while (s.next()) {
        };
        return s.summary().acoustic_samples;
    };
    auto a = std::async(std::launch::async, worker), b = std::async(std::launch::async, worker);
    check(a.get() == b.get(), "Concurrent sample mismatch");
    o.metrics->max_values = 1;
    rejects([&] {
        Simulation s(c, o);
        s.next();
    });
    o.metrics->max_values = 32000000;
    o.metrics->retarded_time = true;
    o.propagation = aeroacoustics::PropagationOptions{};
    o.propagation->ground = aeroacoustics::GroundModel::rigid;
    rejects([&] { Simulation s(c, o); });
}
} // namespace
int main(int argc, char **argv) {
    try {
        check(argc == 4, "metrics_probe CASE SURFACES OUTPUT");
        Case c(argv[1]);
        analytic();
        integration(std::filesystem::path(argv[3]) / "analytic-bins");
        moving_absorption(std::filesystem::path(argv[3]) / "moving-absorption");
        surfaces(argv[2], c);
        simulation(c, argv[2], argv[3]);
        std::cout << "E3/E5/E7 analytic, surface, checkpoint and concurrent checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
