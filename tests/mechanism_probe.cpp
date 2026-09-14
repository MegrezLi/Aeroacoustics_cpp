#include "aeroacoustics.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace aeroacoustics;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
int main() {
    try {
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
        check(tno[index(Mechanism::trailing_separation)] == bpm[index(Mechanism::trailing_separation)],
              "TNO must retain BPM separation");
        const auto direct =
            tblte_tno(s.speed, s.trailing.theta, s.trailing.phi, s.span, s.trailing.distance, s.bl, p);
        check(tno[index(Mechanism::trailing_pressure)] == direct.first &&
                  tno[index(Mechanism::trailing_suction)] == direct.second,
              "TNO must replace, not add BPM pressure/suction");
        std::cout << configurations << " mechanism switch combinations and TNO replacement passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
