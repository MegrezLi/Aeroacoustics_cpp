// Brooks, Pope & Marcolini (1989), NASA-RP-1218, Appendix D / Figure 45(a).
// Reference totals are model outputs from Table D2, NOT experimental data.
#include "aeroacoustics.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        using namespace aeroacoustics;
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout << "Usage: aeroacoustics_bpm1989 [output.csv] [--check]\n"
                         "--check requires max absolute difference from Table D2 <= 0.01 dB.\n";
            return 0;
        }
        if (argc > 3 || (argc == 3 && std::string(argv[2]) != "--check") ||
            (argc > 1 && std::string(argv[1]).rfind("--", 0) == 0))
            throw std::invalid_argument("Usage: aeroacoustics_bpm1989 [output.csv] [--check]");
        const bool check = argc == 3;
        Parameters p;
        p.freqlist = {100,125,160,200,250,315,400,500,630,800,1000,1250,
                      1600,2000,2500,3150,4000,5000,6300,8000,10000,
                      12500,16000,20000,25000,31500,40000};
        const Spectrum reference = {29.336,32.676,36.042,38.815,41.356,
            43.768,46.057,48.034,49.954,51.849,53.568,55.255,57.106,
            58.817,60.167,60.496,59.455,58.208,56.553,54.821,53.192,
            51.523,49.591,47.731,45.737,43.503,40.987};
        p.spdsound = 340.46;
        p.kinvisc = 1.4529e-5;
        p.x_blmethod = 1;
        p.itrip = 0;
        p.tbltemod = 1;
        p.lammod = 1;
        p.timod = p.tipmod = p.bluntmod = 0;
        p.aweighting = false;
        Section s;
        s.chord = .3048;
        s.span = .4572;
        s.speed = 71.3;
        s.alpha_deg = 1.516; // Effective angle, not tunnel geometric angle.
        s.trailing = {1.22, 90., 90.};
        s.is_tip = false;
        s.stall_deg = 12.5; // Library stall threshold; this low-angle case is below it.
        s.tabulated_boundary_layer = false;
        const auto spectra = section_spectrum(p, s);
        const char* filename = argc > 1 ? argv[1] : "bpm1989_comparison.csv";
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot open output CSV");
        out << "frequency_Hz";
        for (const auto& m : mechanism_registry) out << ',' << m.name << "_dB";
        out << ",total_dB,reference_D2_dB,difference_dB\n" << std::setprecision(12);
        double max_error = 0., sum_squared = 0.;
        for (std::size_t i = 0; i < p.freqlist.size(); ++i) {
            Spectrum active;
            out << p.freqlist[i];
            for (const auto& m : mechanism_registry) {
                const double level = spectra[index(m.id)][i];
                out << ',' << level;
                if (m.contributes_to_total) active.push_back(level);
            }
            const double total = db_sum(active);
            if (!std::isfinite(total)) throw std::runtime_error("Non-finite total SPL");
            const double difference = total - reference[i];
            max_error = std::max(max_error, std::abs(difference));
            sum_squared += difference * difference;
            out << ',' << total << ',' << reference[i] << ',' << difference << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("Cannot finish output CSV");
        std::cout << "Output: " << filename << "\nMax absolute difference: "
                  << max_error << " dB\nRMSE: "
                  << std::sqrt(sum_squared / reference.size()) << " dB\n";
        // Regression tolerance for this published-model case, not an experimental accuracy bound.
        if (check && max_error > .01) {
            std::cerr << "Table D2 comparison exceeds 0.01 dB\n";
            return 2;
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
