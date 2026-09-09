// Compile with AERO_BENCHMARK_LEGACY against a prior library for independent
// snapshot comparisons. The default build exercises the reusable C++ API.
#include "aeroacoustics.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace aeroacoustics;

namespace {
std::vector<Node> nodes(std::size_t count, int step) {
    std::vector<Node> result(count);
    for (std::size_t n = 0; n < count; ++n) {
        auto &v = result[n];
        v.section.chord = .4 + .02 * n;
        v.section.span = .7;
        v.section.speed = 40. + n + .2 * step;
        const double angles[] = {-3., 0., 5., 12.5, 22., 64.};
        v.section.alpha_deg = angles[(n + step) % 6];
        v.section.is_tip = n + 1 == count;
        v.section.ti_section = step % 3 == 0 ? 0. : .1;
        v.section.bl = {{{.0012, .0008}}, {{.012, .008}}, {{.003, .002}}, {{.7, 1.1}}};
        if (n % 4 == 0)
            v.section.bl.cf[0] = 0.;
        v.aero_center = {0., double(n), 80. + n};
        v.inflow = {8., 0., 0.};
    }
    return result;
}
std::vector<Vec3> observers(std::size_t count) {
    std::vector<Vec3> result;
    for (std::size_t i = 0; i < count; ++i)
        result.push_back({175. * std::cos(.7 * i), 175. * std::sin(.7 * i), 2.});
    return result;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "verify") {
            std::ofstream out(argv[2]);
            out << std::setprecision(17);
            for (int configuration = 0; configuration < 6; ++configuration) {
                Parameters p;
                p.tbltemod = configuration % 3 == 2 ? 2 : 1;
                p.x_blmethod = configuration % 3 == 0 ? 1 : 2;
                p.itrip = configuration % 3;
                p.timod = configuration % 3;
                p.bluntmod = p.lammod = p.tipmod = 1;
                p.aweighting = configuration >= 3;
                // Non-default frequencies also catch fixed-size assumptions.
                p.freqlist = {10., 100., 1000., 4000., 16000.};
#ifndef AERO_BENCHMARK_LEGACY
                AcousticWorkspace workspace(p);
#endif
                for (int step = 0; step < 6; ++step) {
                    auto input = nodes(step % 2 ? 3 : 7, step);
                    auto receivers = observers(step % 2 ? 2 : 4);
                    if (step == 4)
                        receivers.back() = {1e12, 1e12, 2.};
#ifdef AERO_BENCHMARK_LEGACY
                    const auto snapshot = snapshot_spectrum(p, input, receivers);
#else
                    const auto &snapshot = workspace.evaluate(input, receivers);
                    // Compare reused storage against independent single-section
                    // calls after changing node counts, geometry and flow state.
                    for (std::size_t o = 0; o < receivers.size(); ++o)
                        for (std::size_t n = 0; n < input.size(); ++n) {
                            const auto &node = input[n];
                            const auto geometry = observe(receivers[o], node.aero_center, node.global_to_local,
                                                          node.section.chord, node.airfoil_reference);
                            auto section = node.section;
                            section.leading = geometry.first;
                            section.trailing = geometry.second;
                            const auto expected = section_spectrum(p, section);
                            for (std::size_t m = 0; m < expected.size(); ++m)
                                for (std::size_t f = 0; f < expected[m].size(); ++f) {
                                    const double a = expected[m][f], b = snapshot[o][n][m][f];
                                    if (std::isnan(a) || std::isnan(b) ||
                                        (a != b && (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b)>1e-10)))
                                        throw std::runtime_error("Reused workspace differs from independent section");
                                }
                        }
#endif
                    for (const auto &observer : snapshot)
                        for (const auto &node : observer)
                            for (const auto &mechanism : node)
                                for (double value : mechanism)
                                    out << value << '\n';
                }
            }
            return out ? 0 : 1;
        }
        if (argc != 5)
            throw std::runtime_error(
                "Usage: probe verify OUTPUT | probe benchmark bpm|tno OBSERVERS REPEATS");
        Parameters p;
        p.tbltemod = std::string(argv[2]) == "tno" ? 2 : 1;
        p.bluntmod = 1;
        p.timod = 2;
        auto input = nodes(30, 1);
        auto receivers = observers(std::stoul(argv[3]));
        const int repeats = std::stoi(argv[4]);
        double checksum = 0.;
#ifndef AERO_BENCHMARK_LEGACY
        AcousticWorkspace workspace(p);
#endif
        const auto evaluate = [&] {
#ifdef AERO_BENCHMARK_LEGACY
            const auto snapshot = snapshot_spectrum(p, input, receivers);
#else
            const auto &snapshot = workspace.evaluate(input, receivers);
#endif
            for (const auto &observer : snapshot)
                for (const auto &node : observer)
                    for (const auto &mechanism : node)
                        for (double value : mechanism)
                            if (std::isfinite(value))
                                checksum += value;
        };
        evaluate(); // Warm up allocations and math-library initialization.
        const auto started = std::chrono::steady_clock::now();
        for (int i = 0; i < repeats; ++i)
            evaluate();
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        std::cout << std::setprecision(17) << "{\"seconds\":" << seconds << ",\"checksum\":" << checksum
                  << ",\"repeats\":" << repeats << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
