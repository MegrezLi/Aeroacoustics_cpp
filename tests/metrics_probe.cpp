#include "turbine/simulation.hpp"
#include <algorithm>
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
    EngineeringMetrics free(o, p, {{10, 0, 0}}, 1, {}, {}),
        absorbed(o, p, {{10, 0, 0}}, 1, {}, prop);
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
        double t = plain[2], factor = p.spdsound / (p.spdsound - 1),
               te = (t - 10 / p.spdsound) * factor;
        near(plain[3], 1000 * factor);
        near(plain[4] - loss[4], prop.atmosphere.absorption_db_per_m(plain[3]) * (10 - te));
    }
    // A failed history cannot be exported as a complete result.
    rejects([&] { free.append(0, 8, 1, nodes, 0, 0, spectrum); });
    rejects([&] { free.write(out / "free"); });
}
void surfaces(const std::filesystem::path &file, const Case &c) {
    auto set = SurfaceSet::read(file);
    SurfaceSet copy = *set;
    for (std::size_t i = 0; i < set->data().size(); ++i) {
        const auto id = set->data()[i].airfoil;
        check(set->find(id) == &set->data()[i], "Surface index changed input order");
        check(copy.find(id) == &copy.data()[i], "Surface copy retained original storage");
    }
    check(set->find(std::numeric_limits<std::size_t>::max()) == nullptr, "Unknown surface ID");
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
void simulation(const Case &c, const std::filesystem::path &surface,
                const std::filesystem::path &out) {
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
    check(contents(out / "first/receiver_history.csv") ==
              contents(out / "restored/receiver_history.csv"),
          "Metric checkpoint mismatch");
    check(contents(out / "first/receiver_history.csv") ==
              contents(out / "reset/receiver_history.csv"),
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
namespace regression {
using namespace aeroacoustics;
void relative_near(double a, double b) {
    check(std::isfinite(a) && std::isfinite(b) &&
              std::abs(a - b) <= 1e-10 * std::max({1., std::abs(a), std::abs(b)}),
          "Numeric mismatch");
}
void percentiles() {
    for (double peak : {1e-50, 100., 1e100}) {
        auto s = turbine::level_statistics({0, 9, 10}, {0, 0, peak});
        check(s.l50_db == -INFINITY && s.l95_db == -INFINITY, "Partial silence percentiles");
        relative_near(s.l5_db, 10 * std::log10(.5 * peak));
        relative_near(s.leq_db, 10 * std::log10(.05 * peak));
    }
    auto boundary = turbine::level_statistics({0, .5, 10}, {0, 0, 100});
    check(boundary.l95_db == -INFINITY, "Silence mass exactly at quantile");
    auto below = turbine::level_statistics({0, .4, 10}, {0, 0, 100});
    relative_near(below.l95_db, 10 * std::log10(100. * .1 / 9.6));
    auto constant = turbine::level_statistics({0, 1, 2}, {4, 4, 4});
    relative_near(constant.l95_db, 10 * std::log10(4.));
    auto plateau = turbine::level_statistics({0, 1, 9, 10}, {0, 4, 4, 8});
    relative_near(plateau.l50_db, 10 * std::log10(4.));
    auto silence = turbine::level_statistics({0, 1}, {0, 0});
    check(silence.leq_db == -INFINITY && silence.l5_db == -INFINITY &&
              silence.l50_db == -INFINITY && silence.l95_db == -INFINITY,
          "All-silent statistics");
}
std::vector<std::vector<double>> rows(const std::filesystem::path &file) {
    std::ifstream in(file);
    check(bool(in), "Missing CSV");
    std::string line, word;
    std::getline(in, line);
    std::vector<std::vector<double>> result;
    while (std::getline(in, line)) {
        std::istringstream stream(line);
        std::vector<double> row;
        while (std::getline(stream, word, ','))
            if (!word.empty())
                row.push_back(std::stod(word));
        result.push_back(row);
    }
    return result;
}
double integral(const std::vector<double> &t, const std::vector<double> &p, double begin,
                double end) {
    double result = 0;
    for (std::size_t i = 1; i < t.size(); ++i) {
        double a = std::max(begin, t[i - 1]), b = std::min(end, t[i]);
        if (a >= b)
            continue;
        auto at = [&](double x) {
            return p[i - 1] + (p[i] - p[i - 1]) * (x - t[i - 1]) / (t[i] - t[i - 1]);
        };
        result += (b - a) * (at(a) + at(b)) / 2;
    }
    return result;
}
void metrics(const std::filesystem::path &output) {
    for (bool moving : {false, true}) {
        std::vector<double> reference;
        std::vector<std::vector<double>> reference_bins;
        for (double dt : {.1, .2, .07, .125}) {
            turbine::MetricsOptions options;
            options.start = moving ? .12 : 0;
            options.end = moving ? .36 : .2;
            options.receiver_dt = dt;
            options.am_window = 14;
            options.retarded_time = moving;
            options.background_laeq_db = 0;
            Parameters p;
            p.freqlist = {1000};
            p.aweighting = true;
            const std::size_t n = moving ? 2 : 1;
            turbine::EngineeringMetrics m(options, p, {{100, 0, 0}}, n, {100, 0, 0}, {});
            std::vector<std::vector<Node>> nodes(1, std::vector<Node>(n));
            Snapshot snapshot(1, std::vector<Mechanisms>(n));
            for (auto &node : snapshot[0])
                for (auto &source : node)
                    source = {-INFINITY};
            std::vector<std::vector<double>> times(n), powers(n);
            for (int i = 0; i < (moving ? 5 : 3); ++i) {
                double t = i * .1;
                for (std::size_t j = 0; j < n; ++j) {
                    double energy = j == 0 ? (i == 1 ? 100 : 1) : (i == 2 ? 50 : 2);
                    double delay = moving ? (j == 0 ? .05 + .1 * t : .1 + .2 * t) : 0;
                    nodes[0][j].airfoil_reference = {1, 0};
                    nodes[0][j].aero_center = {100 - p.spdsound * delay, 0, 0};
                    snapshot[0][j][1][0] = 10 * std::log10(energy);
                    times[j].push_back(t + delay);
                    powers[j].push_back(energy);
                }
                m.append(t, i % 2 ? 8.5 : 7.5, 1, nodes, 0, 0, snapshot);
            }
            auto folder = output / ((moving ? "moving-" : "static-") + std::to_string(dt));
            std::filesystem::create_directories(folder);
            m.write(folder);
            auto summary = rows(folder / "receiver_map.csv").at(0);
            double energy = 0;
            for (std::size_t j = 0; j < n; ++j)
                energy += integral(times[j], powers[j], options.start, options.end);
            const double duration = options.end - options.start;
            relative_near(summary[7], 10 * std::log10(energy / duration));
            relative_near(summary[8], 10 * std::log10(energy / duration + 1));
            auto bins = rows(folder / "wind_bins.csv");
            double bin_time = 0, bin_energy = 0;
            for (const auto &bin : bins) {
                bin_time += bin[3];
                bin_energy += bin[3] * std::pow(10., bin[4] / 10);
            }
            relative_near(bin_time, duration);
            relative_near(bin_energy, energy);
            if (reference.empty()) {
                reference = summary;
                reference_bins = bins;
            } else {
                for (std::size_t k = 4; k < summary.size(); ++k)
                    relative_near(summary[k], reference[k]);
                check(bins.size() == reference_bins.size(),
                      "Wind bin count changed with ReceiverDT");
                for (std::size_t j = 0; j < bins.size(); ++j)
                    for (std::size_t k = 0; k < bins[j].size(); ++k)
                        relative_near(bins[j][k], reference_bins[j][k]);
            }
        }
    }
}
} // namespace regression
} // namespace
int main(int argc, char **argv) {
    try {
        check(argc == 4, "metrics_probe CASE SURFACES OUTPUT");
        Case c(argv[1]);
        analytic();
        regression::percentiles();
        regression::metrics(std::filesystem::path(argv[3]) / "energy-grid");
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
