#include "turbine/solver.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
namespace {
double decibels(double p) { return p == 0 ? 0 : 10 * std::log10(p); }
} // namespace
int main(int argc, char **argv) {
    try {
        const auto started = std::chrono::steady_clock::now();
        if (argc < 3 || argc > 4)
            throw std::runtime_error("Usage: aeroacoustics_turbine CASE.fst OUTPUT_DIRECTORY [duration]");
        turbine::Case c(argv[1]);
        if (argc == 4)
            c.duration = std::stod(argv[3]);
        if (!std::isfinite(c.duration) || c.duration < 0)
            throw std::runtime_error("Duration must be finite and nonnegative");
        const std::filesystem::path directory = argv[2];
        std::filesystem::create_directories(directory);
        auto parameters = aeroacoustics::read_aa_input(c.acoustic.path.string()).first;
        parameters.airdens = c.rho;
        parameters.kinvisc = c.nu;
        parameters.spdsound = c.sound_speed;
        const auto observers = aeroacoustics::read_observers(c.acoustic.file("ObserverLocations").string());
        aeroacoustics::Spectrum span;
        for (const auto &s : c.stations)
            span.push_back(s.span);
        const auto elements = aeroacoustics::blade_elements(span, c.acoustic.number("BldPrcnt"));
        const auto first = elements.first;
        aeroacoustics::AcousticDriver acoustic(
            parameters, span, 3, observers, c.acoustic.number("DT_AA", c.dt), c.acoustic.number("AAStart"),
            c.acoustic.number("BldPrcnt"), c.structure.number("TowerHt") + c.structure.number("Twr2Shft"),
            c.acoustic.integer("TICalcMeth"));
        std::vector<std::vector<aeroacoustics::Node>> nodes(3);
        std::vector<std::optional<aeroacoustics::BLTable>> tables(c.stations.size());
        for (int b = 0; b < 3; ++b)
            for (std::size_t j = 0; j < c.stations.size(); ++j) {
                const auto &s = c.stations[j];
                const auto &af = c.airfoils[s.airfoil];
                aeroacoustics::Node n;
                n.section.chord = s.chord;
                n.section.stall_deg = af.input.number("alpha1");
                n.airfoil_reference = af.reference;
                if (parameters.timod == 2) {
                    const auto thickness = aeroacoustics::guidati_thickness(af.coordinates);
                    n.section.thickness_1p = thickness[0];
                    n.section.thickness_10p = thickness[1];
                }
                if (parameters.bluntmod) {
                    turbine::InputFile bl(af.input.file("BL_file"));
                    n.section.te_angle = bl.number("TEAngle");
                    n.section.te_thickness = bl.number("TEThick");
                }
                if ((parameters.x_blmethod == 2 || parameters.tbltemod == 2) && b == 0)
                    tables[j] = aeroacoustics::BLTable::read(af.input.file("BL_file").string());
                nodes[b].push_back(n);
            }
        std::array<std::ofstream, 4> outputs;
        const std::string prefix = std::filesystem::path(c.acoustic.value("AAOutFile")).filename().string();
        std::array<std::vector<std::string>, 4> labels;
        const std::array<std::string, 7> mechanism_names{
            {"LBL", "TBL_pressure", "TBL_suction", "TBL_separation", "bluntness", "tip", "inflow"}};
        for (std::size_t o = 0; o < observers.size(); ++o) {
            const auto observer = "Obs" + std::to_string(o + 1);
            labels[0].push_back(observer);
            for (double f : parameters.freqlist) {
                std::ostringstream value;
                value << f;
                labels[1].push_back(observer + "_Freq" + value.str());
                for (const auto &m : mechanism_names)
                    labels[2].push_back(observer + "_Freq" + value.str() + "_" + m);
            }
        }
        for (int b = 1; b <= 3; ++b)
            for (std::size_t j = 0; j < c.stations.size(); ++j)
                for (std::size_t o = 0; o < observers.size(); ++o)
                    labels[3].push_back("Blade" + std::to_string(b) + "_Node" + std::to_string(j + 1) +
                                        "_Obs" + std::to_string(o + 1));
        for (int k = 0; k < c.acoustic.integer("NrOutFile"); ++k) {
            outputs[k].open(directory / (prefix + std::to_string(k + 1) + ".out"));
            if (!outputs[k])
                throw std::runtime_error("Cannot create acoustic output");
            outputs[k] << "C++ turbine aeroacoustics; reference sound pressure 20 uPa\nTime";
            for (const auto &label : labels[k])
                outputs[k] << '\t' << label;
            outputs[k] << "\n(s)";
            for (std::size_t j = 0; j < labels[k].size(); ++j)
                outputs[k] << "\t(dB)";
            outputs[k] << '\n' << std::setprecision(12);
        }
        std::ofstream dynamics(directory / "dynamics.csv");
        dynamics << std::setprecision(17) << "time";
        for (int b = 1; b <= 3; ++b)
            dynamics << ",q" << b << "_flap1,q" << b << "_flap2,q" << b << "_edge,qd" << b << "_flap1,qd" << b
                     << "_flap2,qd" << b << "_edge";
        dynamics << '\n';
        turbine::Solver solver(c);
        std::size_t snapshots = 0;
        for (;;) {
            dynamics << solver.time;
            for (const auto &s : solver.state)
                for (const auto &v : {s.q, s.qd})
                    for (double x : v)
                        dynamics << ',' << x;
            dynamics << '\n';
            for (int b = 0; b < 3; ++b)
                for (std::size_t j = 0; j < c.stations.size(); ++j) {
                    const auto &a = solver.aerodynamic.blades[b][j];
                    auto &node = nodes[b][j];
                    node.aero_center = a.motion.position;
                    node.inflow = a.wind;
                    node.section.speed = a.speed;
                    node.section.alpha_deg = a.alpha / turbine::deg;
                    for (int i = 0; i < 3; ++i)
                        for (int k = 0; k < 3; ++k)
                            node.global_to_local[3 * i + k] = a.motion.orientation[i][k];
                    if (tables[j])
                        node.section.bl = tables[j]->interpolate(
                            node.section.alpha_deg, a.speed * node.section.chord / c.nu, node.section.chord);
                }
            auto snapshot = acoustic.step(solver.time, nodes);
            if (snapshot) {
                ++snapshots;
                const std::size_t nf = parameters.freqlist.size(), no = observers.size();
                std::vector<double> total(no), spectra(no * nf), mechanisms(no * nf * 7),
                    nodal(3 * c.stations.size() * no);
                for (std::size_t o = 0; o < no; ++o)
                    for (std::size_t n = 0; n < (*snapshot)[o].size(); ++n) {
                        double node_power = 0;
                        for (std::size_t m = 0; m < 7; ++m)
                            for (std::size_t f = 0; f < nf; ++f) {
                                const double db = (*snapshot)[o][n][m][f],
                                             power = std::isfinite(db) ? std::pow(10., db / 10) : 0;
                                total[o] += power;
                                spectra[o * nf + f] += power;
                                mechanisms[(o * nf + f) * 7 + m] += power;
                                node_power += power;
                            }
                        const auto b = n / (c.stations.size() - first),
                                   j = n % (c.stations.size() - first) + first;
                        nodal[(b * c.stations.size() + j) * no + o] = node_power;
                    }
                const std::array<std::vector<double>, 4> values{{total, spectra, mechanisms, nodal}};
                for (int k = 0; k < c.acoustic.integer("NrOutFile"); ++k) {
                    outputs[k] << solver.time;
                    for (double p : values[k])
                        outputs[k] << '\t' << decibels(p);
                    outputs[k] << '\n';
                }
            }
            if (solver.time + c.dt / 2 >= c.duration)
                break;
            solver.step();
        }
        std::ofstream metadata(directory / "run.json");
        metadata << std::setprecision(17) << "{\n  \"solver\": \"standalone C++\",\n  \"dt\": " << c.dt
                 << ",\n  \"duration\": " << solver.time << ",\n  \"steps\": " << solver.step_number
                 << ",\n  \"acoustic_samples\": " << snapshots << ",\n  \"elapsed_seconds\": "
                 << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
                 << "\n}\n";
        std::cout << "C++ turbine completed " << solver.time << " s, " << solver.step_number << " steps, "
                  << snapshots << " acoustic samples\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
