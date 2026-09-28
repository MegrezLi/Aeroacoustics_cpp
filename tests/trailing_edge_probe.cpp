#include "trailing_edge.hpp"
#include "turbine/simulation.hpp"
#include <algorithm>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
using namespace aeroacoustics;
namespace {
void check(bool ok, const char *s) {
    if (!ok)
        throw std::runtime_error(s);
}
void close(double a, double b, double r = 1e-10) {
    check(std::abs(a - b) <= r * std::max({1e-20, std::abs(a), std::abs(b)}),
          "Trailing-edge numeric mismatch");
}
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, "Expected rejection");
}
void api_checks() {
    Parameters p;
    p.timod = 0;
    p.tbltemod = 3;
    p.x_blmethod = 2;
    p.freqlist = {100, 1000, 10000};
    Section s;
    s.speed = 40;
    s.chord = .2;
    s.span = .5;
    s.alpha_deg = 0;
    s.bl = {{{.001, .001}}, {{.007, .007}}, {{.002, .002}}, {{1., 1.}}};
    for (double q : {1e-8, .01, 1., 10., 1e8}) {
        close(howe_chase_shape(q, 0, .04, .007), q * q / std::pow(q * q + 1.33 * 1.33, 2));
        close(howe_chase_shape(q, 1e-12, .04, .007), howe_chase_shape(q, 0, .04, .007), 1e-8);
    }
    const double q = 1e8;
    close(howe_chase_shape(q, .02, .04, .007) / howe_chase_shape(q, 0, .04, .007), 1. / 5, 1e-7);
    const auto straight = section_spectrum(p, s);
    const auto bands = FrequencyBands::openfast_reference(p.freqlist);
    const double pi = std::acos(-1.), chi = 1.33, scale = 2 * pi * .007 / (.7 * s.speed);
    auto primitive = [&](double x) { return std::atan(x / chi) / (2 * chi) - x / (2 * (x * x + chi * chi)); };
    const double coeff = 4 * .1553 * std::pow(p.airdens * .03 * s.speed, 2) * s.span / p.spdsound *
                         std::pow(.007 / s.trailing.distance, 2) * .5;
    for (std::size_t i = 0; i < bands.size(); ++i) {
        const auto &b = bands.values()[i];
        const double power = coeff / scale * (primitive(b.upper_hz * scale) - primitive(b.lower_hz * scale));
        close(std::pow(10., straight[index(Mechanism::trailing_suction)][i] / 10) * 4e-10, power, 1e-9);
        check(straight[index(Mechanism::trailing_separation)][i] == -INFINITY,
              "Howe must not retain BPM separation");
    }
    p.howe.half_height = .02;
    const auto serrated = section_spectrum(p, s);
    check(serrated[1] != straight[1], "Serrations have no effect");
    close(howe_chase_psd(1000, 40, .007, .5, {2, 90, 90}, p) * 4,
          howe_chase_psd(1000, 40, .007, .5, {1, 90, 90}, p));
    close(howe_chase_psd(1000, 40, .007, 1, {1, 90, 90}, p),
          2 * howe_chase_psd(1000, 40, .007, .5, {1, 90, 90}, p));
    check(howe_chase_psd(1000, 40, .007, .5, {1, 0, 90}, p) == 0, "Howe null directivity");
    rejects([&] { howe_chase_psd(1000, 150, .007, .5, {1, 90, 90}, p); });
    auto invalid = s;
    invalid.bl.cf[0] = 0;
    rejects([&] { section_spectrum(p, invalid); });
    invalid = s;
    invalid.alpha_deg = 30;
    rejects([&] { section_spectrum(p, invalid); });
    rejects([&] { howe_chase_shape(1, .02, 0, .007); });
    // Same geometry and data must give the same answers through the workspace and scalar APIs.
    Node node;
    node.section = s;
    const Vec3 observer{1, 2, 3};
    const auto geometry = observe(observer, node.aero_center, node.global_to_local, s.chord);
    s.leading = geometry.first;
    s.trailing = geometry.second;
    for (int model : {2, 3})
        for (auto mode : {TnoEdgeVelocity::reference, TnoEdgeVelocity::input}) {
            p.tbltemod = model;
            p.tno_edge_velocity = mode;
            s.bl.edge_velocity_ratio = {1.15, .85};
            node.section = s;
            const auto scalar = section_spectrum(p, s);
            AcousticWorkspace workspace(p);
            const auto result = workspace.evaluate({node}, {observer})[0][0];
            for (std::size_t m = 0; m < mechanism_count; ++m)
                for (std::size_t f = 0; f < p.freqlist.size(); ++f)
                    if (std::isfinite(scalar[m][f]))
                        close(result[m][f], scalar[m][f]);
                    else
                        check(result[m][f] == scalar[m][f], "Inactive channel mismatch");
            auto worker = [=] { return section_spectrum(p, s); };
            auto task = std::async(std::launch::async, worker);
            check(task.get() == scalar, "Concurrent trailing-edge result differs");
        }
    p.tbltemod = 2;
    p.tno_edge_velocity = TnoEdgeVelocity::reference;
    const auto reference = section_spectrum(p, s);
    s.bl.edge_velocity_ratio = {1, 1};
    check(reference == section_spectrum(p, s), "Reference TNO must use unit edge ratios");
    p.tno_edge_velocity = TnoEdgeVelocity::input;
    check(reference == section_spectrum(p, s), "Unit input TNO differs from reference");
    s.bl.edge_velocity_ratio = {1.15, -.85};
    auto varied = section_spectrum(p, s);
    s.bl.edge_velocity_ratio[1] = .85;
    check(varied == section_spectrum(p, s), "TNO input sign convention differs");
    check(varied[1] != reference[1] && varied[2] != reference[2], "TNO ignored input ratios");
    auto direct =
        configured_tblte_tno(s.speed, s.trailing.theta, s.trailing.phi, s.span, s.trailing.distance, s.bl, p);
    check(varied[1] == direct.first && varied[2] == direct.second, "Configured TNO scalar differs");
    p.tno_edge_velocity = TnoEdgeVelocity::reference;
    const auto raw =
        tblte_tno(s.speed, s.trailing.theta, s.trailing.phi, s.span, s.trailing.distance, s.bl, p);
    check(raw == direct, "Legacy raw TNO must retain input semantics");
    p.tno_edge_velocity = TnoEdgeVelocity::input;
    s.bl.edge_velocity_ratio[0] = 0;
    rejects([&] { section_spectrum(p, s); });
    s.bl.edge_velocity_ratio[0] = NAN;
    rejects([&] { section_spectrum(p, s); });
}
void write_reference_rows(const std::filesystem::path &folder) {
    std::filesystem::create_directory(folder);
    std::ofstream shapes(folder / "howe.csv"), tno(folder / "tno.csv");
    shapes << std::setprecision(17) << "q,h,lambda,delta,shape\n";
    for (double q : {.01, .1, 1., 10., 100., 1000.})
        for (double h : {0., .001, .01, .02, .1})
            for (double lambda : {.005, .02, .1})
                shapes << q << ',' << h << ',' << lambda << ",0.007," << howe_chase_shape(q, h, lambda, .007)
                       << '\n';
    tno << std::setprecision(17) << "mode,suction_ratio,pressure_ratio,frequency,pressure,suction\n";
    Parameters p;
    p.tbltemod = 2;
    BoundaryLayer bl{{{.0012, .0008}}, {{.0110586, .00746583}}, {{.000378576, .00198438}}, {{1., 1.}}};
    for (auto ratio : std::vector<std::array<double, 2>>{{1, 1}, {1.1, .8}, {-1.1, -.8}})
        for (auto mode : {TnoEdgeVelocity::reference, TnoEdgeVelocity::input}) {
            p.tno_edge_velocity = mode;
            bl.edge_velocity_ratio = ratio;
            auto pair = configured_tblte_tno(63.92, 90, 90, .509, 1.22, bl, p);
            for (std::size_t i = 0; i < p.freqlist.size(); ++i)
                tno << tno_edge_velocity_name(mode) << ',' << ratio[0] << ',' << ratio[1] << ','
                    << p.freqlist[i] << ',' << pair.first[i] << ',' << pair.second[i] << '\n';
        }
    shapes.close();
    tno.close();
    check(bool(shapes) && bool(tno), "Reference CSV write failed");
}
} // namespace
int main(int argc, char **argv) {
    try {
        api_checks();
        if (argc > 1)
            write_reference_rows(argv[1]);
        if (argc > 3) {
            turbine::Case c(argv[2]);
            turbine::RunOptions o;
            o.duration = .2;
            o.surfaces = turbine::SurfaceSet::read(std::filesystem::path(argv[3]) / "surfaces.dat");
            o.trailing_edge =
                TrailingEdgeOptions::read((std::filesystem::path(argv[3]) / "howe-serrated.dat").string());
            turbine::Simulation simulation(c, o);
            simulation.next();
            const auto saved = simulation.checkpoint();
            auto sample = [&] {
                while (auto frame = simulation.next())
                    if (frame->acoustics)
                        return frame->acoustics->power;
                throw std::runtime_error("Missing acoustic sample");
            };
            auto original = sample();
            simulation.restore(saved);
            auto restored = sample();
            check(original == restored, "Trailing-edge checkpoint differs");
            simulation.reset();
            simulation.next();
            check(original == sample(), "Trailing-edge reset differs");
        }
        std::cout
            << "E4/E9 spectral limits, SI bands, configuration, workspace, state and concurrency passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
