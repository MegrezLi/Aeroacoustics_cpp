#include "acoustic_levels.hpp"
#include "aeroacoustics.hpp"
#include "aeroacoustics_c.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace aeroacoustics;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
namespace {
void relative_near(double a, double b) {
    check(std::isfinite(a) && std::isfinite(b) &&
              std::abs(a - b) <= 1e-10 * std::max({1., std::abs(a), std::abs(b)}),
          "Numeric mismatch");
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
void sums() {
    for (auto invalid : {NAN, INFINITY})
        for (Spectrum s : {Spectrum{-INFINITY, invalid}, Spectrum{invalid, -INFINITY},
                           Spectrum{20, invalid, 30}})
            rejects([&] { db_sum(s); });
    check(db_sum({}) == -INFINITY && db_sum({-INFINITY, -INFINITY}) == -INFINITY, "Silence sum");
    relative_near(db_sum({0, 0, -INFINITY}), 10 * std::log10(2.));
    relative_near(db_sum({4000, 4000}), 4000 + 10 * std::log10(2.));
}
void abi() {
    const std::array<double, 26> valid{{3,     .2,   40, 90,   90,   .5,   10,  12,  .007,
                                        .001,  .001, .1, .001, 10,   .02,  .12, 340, 1.48e-5,
                                        1.225, 40,   1,  .002, .002, .007, 1,   1}};
    const std::array<int, 3> flags{{2, 1, 1}};
    const std::array<double, 2> frequencies{{100, 1000}};
    auto run = [&](int op, std::array<double, 26> x, std::array<int, 3> modes,
                   std::array<double, 2> freq, bool fail) {
        std::array<double, 6> out;
        out.fill(1234567.);
        auto before = aeroacoustics_call_count();
        int status = aeroacoustics_kernel(op, 2, freq.data(), x.data(), modes.data(), out.data());
        check((status != 0) == fail, "Wrong C ABI status");
        if (fail) {
            check(std::string(aeroacoustics_last_error()).size() > 0, "Missing C ABI diagnostic");
            check(aeroacoustics_call_count() == before, "Failed call incremented counter");
            for (double v : out)
                check(v == 1234567., "Failure partially wrote output");
        } else {
            check(std::string(aeroacoustics_last_error()).empty(), "Success retained stale error");
            check(aeroacoustics_call_count() == before + 1, "Successful call counter");
            int channels = op == 2 ? 3 : op == 7 ? 2 : 1;
            for (int i = 0; i < 2 * channels; ++i)
                validate_level(out[i]);
        }
        return out;
    };
    for (int op = 1; op <= 7; ++op) {
        for (double v : {double(NAN), double(INFINITY), -1., 0.}) {
            auto f = frequencies;
            f[1] = v;
            run(op, valid, flags, f, true);
        }
        for (double v : {double(NAN), double(INFINITY), -1., 0.}) {
            auto x = valid;
            x[2] = v;
            run(op, x, flags, frequencies, true);
        }
        if (op != 6) {
            auto x = valid;
            x[6] = 0;
            run(op, x, flags, frequencies, true);
            x = valid;
            x[2] = x[16];
            run(op, x, flags, frequencies, true);
        }
        run(op, valid, flags, frequencies, false);
    }
    for (auto pair : std::vector<std::pair<int, int>>{
             {1, 8}, {2, 9}, {3, 20}, {4, 11}, {5, 12}, {6, 14}, {7, 24}}) {
        auto x = valid;
        x[pair.second] = NAN;
        run(pair.first, x, flags, frequencies, true);
    }
    for (int op : {1, 2, 5}) {
        auto invalid = flags;
        invalid[0] = 0;
        run(op, valid, invalid, frequencies, true);
    }
    auto invalid = flags;
    invalid[2] = 2;
    run(3, valid, invalid, frequencies, true);
    auto x = valid;
    x[13] = 15;
    run(5, x, flags, frequencies, true);
    x = valid;
    x[11] = -1;
    run(4, x, flags, frequencies, true);
    x = valid;
    x[11] = 0;
    run(4, x, flags, frequencies, false);
    x = valid;
    x[20] = 1e308;
    run(3, x, flags, frequencies, true); // finite inputs, invalid computed result
    x = valid;
    for (int i : {0, 1, 7, 9, 10, 11, 12, 13, 14, 15, 19, 20})
        x[i] = NAN;
    auto reference = run(7, x, flags, frequencies, false); // unused fields must not be inspected
    x[24] = -1;
    x[25] = -1;
    check(reference == run(7, x, flags, frequencies, false), "TNO signed ratio compatibility");
    x[24] = 1.15;
    x[25] = .85;
    check(reference != run(7, x, flags, frequencies, false), "TNO input ratios ignored");
    x[21] = 0;
    x[22] = -1;
    x[24] = x[25] = 0;
    x[8] = x[23] = 0;
    run(7, x, flags, frequencies, false); // legacy nonpositive Cf disables a side
    x = valid;
    x[7] = x[9] = x[10] = NAN;
    run(1, x, {{2, -1, -1}}, frequencies, false); // tabulated laminar kernel uses only d99P
    x.fill(NAN);
    x[1] = .2;
    x[2] = 40;
    x[14] = .02;
    x[15] = .12;
    run(6, x, {{-1, -1, -1}}, frequencies, false);
    std::array<double, 6> out{};
    check(aeroacoustics_kernel(4, 0, frequencies.data(), valid.data(), flags.data(), out.data()) !=
              0,
          "Zero n");
    check(aeroacoustics_kernel(4, 2, nullptr, valid.data(), flags.data(), out.data()) != 0,
          "Null frequency");
    run(0, valid, flags, frequencies, true);
    run(4, valid, flags, frequencies, false);
}
} // namespace
int main() {
    try {
        sums();
        abi();
        const char *legacy[] = {"LBL",       "TBL_pressure", "TBL_suction", "TBL_separation",
                                "bluntness", "tip",          "inflow"};
        check(mechanism_count == 7, "Legacy channel count changed");
        for (std::size_t i = 0; i < mechanism_count; ++i)
            check(std::string(mechanism_registry[i].name) == legacy[i] &&
                      index(mechanism_registry[i].id) == i,
                  "Legacy mechanism mapping changed");
        Section s;
        s.bl = {{{.0012, .0008}}, {{.012, .008}}, {{.003, .002}}, {{1., 1.}}};
        Parameters p;
        p.freqlist = {100, 1000};
        p.x_blmethod = 2;
        std::size_t configurations = 0;
        for (int trailing = 0; trailing < 3; ++trailing)
            for (int inflow = 0; inflow < 3; ++inflow)
                for (int trip = 0; trip < 3; ++trip)
                    for (int laminar = 0; laminar < 2; ++laminar)
                        for (int tip = 0; tip < 2; ++tip)
                            for (int blunt = 0; blunt < 2; ++blunt) {
                                p.tbltemod = trailing;
                                p.timod = inflow;
                                p.itrip = trip;
                                p.lammod = laminar;
                                p.tipmod = tip;
                                p.bluntmod = blunt;
                                const auto out = section_spectrum(p, s);
                                const bool enabled[] = {bool(laminar && trip == 0),
                                                        bool(trailing),
                                                        bool(trailing),
                                                        bool(trailing),
                                                        bool(blunt),
                                                        bool(tip),
                                                        bool(inflow)};
                                for (std::size_t m = 0; m < mechanism_count; ++m)
                                    for (double value : out[m])
                                        check(enabled[m] ? !std::isnan(value) && value != INFINITY
                                                         : value == -INFINITY,
                                              "Mechanism activation mismatch");
                                ++configurations;
                            }
        p.tbltemod = 2;
        p.timod = p.tipmod = p.bluntmod = p.lammod = 0;
        auto tno = section_spectrum(p, s);
        p.tbltemod = 1;
        auto bpm = section_spectrum(p, s);
        check(tno[index(Mechanism::trailing_separation)] ==
                  bpm[index(Mechanism::trailing_separation)],
              "TNO must retain BPM separation");
        const auto direct = tblte_tno(s.speed, s.trailing.theta, s.trailing.phi, s.span,
                                      s.trailing.distance, s.bl, p);
        check(tno[index(Mechanism::trailing_pressure)] == direct.first &&
                  tno[index(Mechanism::trailing_suction)] == direct.second,
              "TNO must replace, not add BPM pressure/suction");
        std::cout << configurations
                  << " mechanism switch combinations and TNO replacement passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
