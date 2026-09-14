#include "lookup_diagnostics.hpp"
#include "turbine/solver.hpp"
#include <cstdlib>
#include <iostream>
#include <new>

namespace {
bool counting = false;
std::size_t allocations = 0;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
template <class F> void rejects(F f) {
    bool rejected = false;
    try {
        f();
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "Expected invalid input rejection");
}
} // namespace
void *operator new(std::size_t size) {
    if (counting)
        ++allocations;
    if (void *p = std::malloc(size ? size : 1))
        return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }

int main(int argc, char **argv) {
    try {
        check(argc == 2, "Usage: turbine_workspace_probe CASE.fst");
        turbine::Case c(argv[1]);
        turbine::Solver solver(c);
        solver.step(); // Initial sizing of every Newton and mapping buffer.
        counting = true;
        for (int i = 0; i < 50; ++i)
            solver.step();
        counting = false;
        check(allocations == 0, "Solver step allocated after warm-up");

        // Known rigid translation and constant line load, with dirty output buffers.
        const turbine::Matrix3 identity{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
        turbine::MotionMap motion_map({{{0, 0, 0}, identity}, {{2, 0, 0}, identity}},
                                      {{{1, 0, 0}, identity}});
        std::vector<turbine::Motion> source_motion(2), mapped_motion(7);
        for (int step = 0; step < 2; ++step) {
            for (int j = 0; j < 2; ++j) {
                source_motion[j].position = {double(2 * j + step), 0, 0};
                source_motion[j].orientation = identity;
                source_motion[j].velocity = {double(step), 0, 0};
            }
            motion_map.transfer_into(source_motion, mapped_motion);
            check(mapped_motion.size() == 1 && mapped_motion[0].position[0] == step + 1 &&
                      mapped_motion[0].velocity[0] == step,
                  "Motion buffer retained old state");
        }
        rejects([&] { motion_map.transfer_into(source_motion, source_motion); });
        const std::vector<turbine::Vec3> source_positions{{0, 0, 0}, {2, 0, 0}}, destination{{1, 0, 0}};
        turbine::LoadMap load_map(source_positions, destination);
        std::vector<turbine::PointLoad> distributed(2), mapped_loads(7);
        for (auto &load : distributed)
            load.force = {2, 0, 0};
        load_map.transfer_into(distributed, source_positions, destination, mapped_loads);
        check(mapped_loads.size() == 1 && mapped_loads[0].force[0] == 4 &&
                  turbine::norm(mapped_loads[0].moment) == 0,
              "Constant line load not conserved");
        distributed.assign(2, {});
        load_map.transfer_into(distributed, source_positions, destination, mapped_loads);
        check(turbine::norm(mapped_loads[0].force) == 0 && turbine::norm(mapped_loads[0].moment) == 0,
              "Load buffer retained old force/moment");
        rejects([&] { load_map.transfer_into(distributed, source_positions, destination, distributed); });

        using namespace aeroacoustics;
        BLTable table;
        table.aoa = {-10, 10};
        table.reynolds = {1e6, 3e6};
        for (double re : table.reynolds)
            for (double a : table.aoa) {
                std::array<double, 8> row{};
                for (int k = 0; k < 8; ++k)
                    row[k] = k + 2 + a / 10 + re / 1e6;
                table.values.push_back(row);
            }
        auto prepared = table.prepare();
        const auto midpoint = prepared.interpolate(0, 2e6, 2);
        check(midpoint.dstar[0] == 12 && midpoint.dstar[1] == 14 && midpoint.d99[0] == 16 &&
                  midpoint.cf[0] == 10 && midpoint.edge_velocity_ratio[0] == 4,
              "Bilinear interpolation/units");
        check(prepared.interpolate(-20, 0, 1).dstar[0] == 4, "Lower clamp");
        check(prepared.interpolate(20, 4e6, 1).dstar[0] == 8, "Upper clamp");
        table.aoa[1] = table.aoa[0];
        table.values.clear();
        check(prepared.interpolate(0, 2e6, 2).dstar[0] == 12, "Prepared table aliases mutable input");
        rejects([&] { table.interpolate(0, 2e6, 2); });
        rejects([&] { table.prepare(); });
        rejects([&] { prepared.interpolate(0, 2e6, 0); });
        {
            diagnostics::LookupReport report;
            diagnostics::LookupSession session(report);
            prepared.interpolate(-20, 2e6, 1);
            check(report.calls == 1, "Prepared BL lost lookup diagnostics");
            report.policy = diagnostics::LookupPolicy::error;
            bool rejected = false;
            try {
                prepared.interpolate(0, 4e6, 1);
            } catch (const std::exception &) {
                rejected = true;
            }
            check(rejected, "Prepared BL lost strict lookup policy");
        }

        Parameters p;
        p.freqlist = {100, 1000};
        const Spectrum span{0, 2, 4, 6};
        std::vector<std::vector<Node>> nodes(1, std::vector<Node>(span.size()));
        for (std::size_t j = 0; j < span.size(); ++j) {
            nodes[0][j].aero_center = {0, span[j], 80};
            nodes[0][j].inflow = {8, 0, 0};
        }
        for (int method : {1, 2}) {
            AcousticDriver driver(p, span, 1, {{100, 0, 2}}, .1, .2, 70, 80, method);
            TurbulenceState reference(span, 1, .1, 80, method, p.ti, p.avgv);
            const auto elements = blade_elements(span, 70);
            check(driver.first_node() == elements.first, "Selected node start");
            for (int step = 0; step <= 64; ++step) {
                const double time = step * .00625;
                Spectrum speeds;
                std::vector<Vec3> inflow, leading;
                std::vector<Node> selected;
                for (std::size_t j = 0; j < span.size(); ++j) {
                    auto &node = nodes[0][j];
                    node.section.speed = 40 + step + j;
                    node.inflow[0] = 8 + .1 * std::sin(step);
                    speeds.push_back(node.section.speed);
                    inflow.push_back(node.inflow);
                    auto position = node.aero_center;
                    position[1] -= .25 * node.section.chord;
                    leading.push_back(position);
                    if (j >= elements.first) {
                        auto copy = node;
                        copy.section.ti_section = reference.values[j];
                        copy.section.span = elements.second[j];
                        copy.section.is_tip = j + 1 == span.size();
                        selected.push_back(copy);
                    }
                }
                const bool expected_sample = step >= 32 && step % 16 == 0;
                check(driver.is_sample_time(time) == expected_sample, "Sampling schedule");
                const auto *snapshot = driver.step_view(time, nodes);
                check(bool(snapshot) == expected_sample, "Sampler and step disagree");
                if (snapshot)
                    check(*snapshot == snapshot_spectrum(p, selected, {{100, 0, 2}}),
                          "Acoustic sample must use previous TI");
                reference.update(speeds, inflow, leading);
                check(driver.state.values == reference.values, "TI missed non-sampling update");
            }
        }
        std::cout << "P1-P3 checks passed; allocations in 50 warmed solver steps: " << allocations << '\n';
    } catch (const std::exception &e) {
        counting = false;
        std::cerr << e.what() << '\n';
        return 1;
    }
}
