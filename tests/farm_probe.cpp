#include "turbine/farm.hpp"
#include <future>
#include <iostream>
using namespace turbine;
void check(bool ok) {
    if (!ok)
        throw std::runtime_error("E6 analytic check failed");
}
void close(double a, double b) { check(std::abs(a - b) < 1e-9 * std::max(1., std::abs(b))); }
template <class F> void rejects(F f) {
    bool bad = false;
    try {
        f();
    } catch (const std::exception &) {
        bad = true;
    }
    check(bad);
}
int main(int argc, char **argv) {
    try {
        check(argc == 2);
        Case c(argv[1]);
        TowerInfluence tower;
        tower.provenance = "Synthetic test";
        tower.stations = {{0, 4, 1}, {100, 4, 1}};
        tower.validate();
        close(tower.apply({10, 0, 0}, {-10, 0, 50})[0], 9.6);
        close(tower.apply({10, 0, 0}, {0, 10, 50})[0], 10.4);
        close(tower.apply({0, 10, 0}, {0, -10, 50})[1], 9.6);
        close(tower.apply({10, 0, 0}, {10, 0, 103})[0], 10);
        rejects([&] { tower.apply({10, 0, 0}, {1, 0, 50}); });
        tower.potential = false;
        tower.shadow = true;
        close(tower.apply({10, 0, 0}, {20, 0, 50})[0], 10 * (1 - 1 / std::sqrt(10.)));
        close(tower.apply({10, 0, 0}, {-20, 0, 50})[0], 10);
        close(tower.apply({10, 0, 0}, {2.1, 0, 50})[0], 5);
        close(jensen_deficit(.75, 10, 50, .05, 500, 0), 5 * std::pow(50. / 75, 2));
        close(jensen_deficit(.75, 10, 50, .05, -1, 0), 0);
        close(jensen_deficit(.75, 10, 50, .05, 500, 80), 0);
        rejects([] { jensen_deficit(1, 10, 50, .05, 500, 0); });
        auto history = std::make_shared<WakeHistory>();
        history->hub = {0, 0, 100};
        history->radius = 50;
        history->samples = {{0, 10, .5}, {1, 8, .7}};
        history->validate();
        close(history->at(.5).ct, .6);
        close(history->at(.5).inflow_speed, 9);
        rejects([&] { history->at(2); });
        auto wind = c.wind;
        wind.speed = 10;
        wind.exponent = 0;
        wind.propagation = 0;
        wind.upflow = 0;
        FarmWind field(wind, {0, 0, 0}, .05, {history});
        close(field.at(.5, {500, 0, 100})[0], 10 - jensen_deficit(.6, 9, 50, .05, 500, 0));
        FarmWind shifted(wind, {100, 0, 0}, .05, {history});
        close(shifted.at(.5, {400, 0, 100})[0], field.at(.5, {500, 0, 100})[0]);
        wind.propagation = pi / 2;
        FarmWind rotated(wind, {0, 0, 0}, .05, {history});
        close(rotated.at(.5, {0, -500, 100})[1], -field.at(.5, {500, 0, 100})[0]);
        const auto bands = aeroacoustics::FrequencyBands::openfast_reference({1000});
        StationarySource source;
        source.name = "test";
        source.provenance = "synthetic";
        source.kind = "tones";
        source.frequencies = {1050}; // Off-center line: absorption must not use the 1000 Hz band center.
        source.levels_db = {80};
        source.validate();
        const double a = source.receiver_power({10, 0, 0}, bands, 1, 400, {})[0],
                     b = source.receiver_power({20, 0, 0}, bands, 1, 400, {})[0];
        close(a / b, 4);
        aeroacoustics::PropagationOptions absorption;
        const double attenuated = source.receiver_power({10, 0, 0}, bands, 1, 400, absorption)[0];
        close(10 * std::log10(a / attenuated), 10 * absorption.atmosphere.absorption_db_per_m(1050));
        close(10 * std::log10((a + a) / a), 10 * std::log10(2.));
        source.directivity = 1;
        close(source.receiver_power({10, 0, 0}, bands, 1, 400, {})[0], 2 * a);
        close(source.receiver_power({-10, 0, 0}, bands, 1, 400, {})[0], 0);
        rejects([&] { source.receiver_power({0, 0, 0}, bands, 1, 400, {}); });
        Rotor rotor(c);
        auto output = rotor.evaluate(0, {});
        double length = 0;
        for (std::size_t j = 1; j < c.stations.size(); ++j) {
            const auto &a = c.stations[j - 1];
            const auto &b = c.stations[j];
            length += norm(Vec3{b.curve - a.curve, b.sweep - a.sweep, b.span - a.span});
        }
        for (auto &blade : output.blades)
            for (auto &station : blade)
                station.load.force = {2, 0, 0};
        close(rotor.aerodynamic_thrust(output, {1, 0, 0}), 6 * length);
        RunOptions options;
        options.duration = .1;
        auto influence = std::make_shared<TowerInfluence>(tower);
        influence->potential = true;
        influence->shadow = false;
        options.solver.tower = influence;
        Simulation sim(c, options);
        sim.next();
        auto saved = sim.checkpoint();
        auto next = sim.next();
        const double first = next->state[0].q[0];
        sim.restore(saved);
        close(sim.next()->state[0].q[0], first);
        sim.reset();
        sim.next();
        close(sim.next()->state[0].q[0], first);
        auto worker = [&]() {
            Simulation independent(c, options);
            independent.next();
            return independent.next()->state[0].q[0];
        };
        auto task1 = std::async(std::launch::async, worker), task2 = std::async(std::launch::async, worker);
        close(task1.get(), first);
        close(task2.get(), first);
        std::cout << "E6 tower, Jensen, coordinates, sound-power propagation, thrust, checkpoint and "
                     "concurrent checks "
                     "passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
