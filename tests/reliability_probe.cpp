#include "acoustic_levels.hpp"
#include "checked_output.hpp"
#include "lookup_diagnostics.hpp"
#include "turbine/input.hpp"
#include <iostream>
#include <sstream>
#include <streambuf>

namespace {
void expect(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
template <class F> std::string failure(F f) {
    try {
        f();
    } catch (const std::exception &e) {
        return e.what();
    }
    throw std::runtime_error("Expected failure was not raised");
}
struct BadSync : std::stringbuf {
    int sync() override { return -1; }
};
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 3)
            throw std::runtime_error("Usage: reliability_probe AIRFOIL OUTPUT_DIRECTORY");
        using namespace aeroacoustics;
        const auto infinity = std::numeric_limits<double>::infinity();
        const auto nan = std::numeric_limits<double>::quiet_NaN();
        expect(relative_power(-infinity) == 0, "Silence must be zero energy");
        expect(relative_power(0) == 1, "Physical 0 dB must retain positive energy");
        expect(output_decibels(0) == 0 && output_decibels(1) == 0, "Legacy encoding changed");
        for (double value : {nan, infinity, 4000.})
            failure([&] { relative_power(value); });
        for (double value : {nan, infinity, -1.})
            failure([&] { output_decibels(value); });

        // A NaN generated inside a model must reach callers with model/location context.
        Parameters p;
        p.tbltemod = 0;
        p.timod = 0;
        p.tipmod = 1;
        p.alprat = 1e308;
        Node node;
        node.section.alpha_deg = 3.;
        node.aero_center = {0, 0, 80};
        node.blade_number = 2;
        node.node_number = 9;
        AcousticWorkspace workspace(p);
        const auto error = failure([&] { workspace.evaluate({node}, {{175, 0, 2}}); });
        expect(error.find("observer=1") != std::string::npos && error.find("blade=2") != std::string::npos &&
                   error.find("node=9") != std::string::npos &&
                   error.find("mechanism=tip") != std::string::npos &&
                   error.find("frequency_Hz=") != std::string::npos,
               "Missing acoustic error context");

        diagnostics::LookupReport report;
        {
            diagnostics::LookupSession session(report);
            diagnostics::LookupLocation location({2.5, 2, 9, "test"});
            BLTable table;
            table.aoa = {-1, 1};
            table.reynolds = {1e6, 2e6};
            table.values.resize(4, {{1, 1, .01, .02, .03, .04, .005, .006}});
            const auto endpoint = table.interpolate(1, 2e6, 2);
            expect(report.calls == 0, "Exact endpoints must not count as out of range");
            const auto clamped = table.interpolate(2, 3e6, 2);
            expect(clamped.dstar == endpoint.dstar && report.calls == 2, "BL clamp or counts changed");
            table.interpolate(-2, .5e6, 2);
            expect(report.calls == 4 && report.entries.size() == 2, "BL axes not counted separately");
            for (const auto &entry : report.entries)
                expect(entry.second.calls == 2 && entry.second.first_time == 2.5, "Wrong lookup summary");
            turbine::Airfoil airfoil(argv[1]);
            airfoil.alpha = {-.1, .1};
            airfoil.coefficients = {{1, 2, 3}, {4, 5, 6}};
            expect(airfoil.at(.1).cl == 4 && report.calls == 4, "Airfoil endpoint changed");
            expect(airfoil.at(.2).cl == 4 && report.calls == 5, "Airfoil clamp not recorded");
            report.policy = diagnostics::LookupPolicy::error;
            const auto range_error = failure([&] { airfoil.at(.2); });
            expect(range_error.find("time=2.5") != std::string::npos &&
                       range_error.find("blade=2") != std::string::npos &&
                       range_error.find("axis=alpha_deg") != std::string::npos,
                   "Missing lookup context");
            failure([&] { table.interpolate(0, 3e6, 2); });
            failure([&] { airfoil.at(nan); });
            failure([&] { table.interpolate(0, 1e6, infinity); });
            diagnostics::LookupReport nested;
            {
                diagnostics::LookupSession inner(nested);
                diagnostics::check_lookup("inner", "x", 2, 0, 1);
            }
            expect(diagnostics::active_report == &report && nested.calls == 1,
                   "Nested session did not restore context");
        }
        expect(diagnostics::active_report == nullptr, "Lookup session leaked");

        const std::filesystem::path output(argv[2]);
        std::filesystem::create_directories(output);
        diagnostics::CheckedOutput stream;
        stream.open(output / "checked.txt");
        stream << "ok" << 17;
        stream.finish();
        std::ifstream read(output / "checked.txt");
        std::string value;
        read >> value;
        expect(value == "ok17", "Checked output content changed");
        diagnostics::CheckedOutput bad;
        const auto open_error = failure([&] { bad.open(output); });
        expect(open_error.find(output.string()) != std::string::npos, "Missing failed output path");
        BadSync buffer;
        std::ostream sink(&buffer);
        sink << "buffered data";
        const auto flush_error = failure([&] { diagnostics::flush_checked(sink, "injected-failure"); });
        expect(flush_error.find("injected-failure") != std::string::npos, "Flush failure lost");
#ifndef _WIN32
        if (std::filesystem::exists("/dev/full")) {
            diagnostics::CheckedOutput full;
            full.open("/dev/full");
            failure([&] {
                full << std::string(100000, 'x');
                full.finish();
            });
        }
#endif
        std::cout << "R1-R3 reliability checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
