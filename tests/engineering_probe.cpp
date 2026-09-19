#include "turbine/simulation.hpp"
#include <future>
#include <iostream>
using namespace turbine;
namespace {
void check(bool b, const char *m) {
    if (!b)
        throw std::runtime_error(m);
}
void near(double a, double b, double tol = 1e-10) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tol)
        throw std::runtime_error("Engineering numerical check: " + std::to_string(a) +
                                 " != " + std::to_string(b));
}
template <class F> void rejects(F f) {
    bool caught = false;
    try {
        f();
    } catch (const std::exception &) {
        caught = true;
    }
    check(caught, "Expected invalid engineering input to fail");
}
void wind() {
    std::array<std::vector<double>, 4> a{{{0, 2}, {-1, 3}, {-2, 2}, {10, 100}}};
    std::vector<Vec3> v;
    for (auto t : a[0])
        for (auto x : a[1])
            for (auto y : a[2])
                for (auto z : a[3])
                    v.push_back({8 + t + .1 * x + .2 * y + .01 * z, -t + y, .05 * z});
    GridWind f(a, v);
    auto p = f.at(.7, {.3, .4, 60});
    near(p[0], 9.41);
    near(p[1], -.3);
    near(p[2], 3);
    rejects([&] { f.at(3, {0, 0, 50}); });
    rejects([&] { f.at(1, {0, 0, 101}); });
    GridWind uniform({{{0}, {0}, {0}, {0}}}, {{8, 1, 2}});
    near(uniform.at(100, {20, -9, 1000})[0], 8);
}
void drivetrain(ControlConfig c) {
    c.rotor_inertia = 4;
    c.generator_inertia = .5;
    c.gear_ratio = 2;
    c.shaft_stiffness = 12;
    c.shaft_damping = 0;
    c.optimal_torque_gain = 1e-30;
    c.rated_power = 1e-20;
    c.pitch_kp = c.pitch_ki = 0;
    c.noise_start = 100;
    c.max_step = .0005;
    OperatingController ctrl(c, 1, c.pitch_min, 0);
    for (int i = 0; i < 100; ++i)
        ctrl.advance(.01, 1, 0);
    const double jr = 4, jg = 2, wn = std::sqrt(12 * (1 / jr + 1 / jg)), t = 1,
                 relative = std::sin(wn * t) / (jr * wn);
    const auto &s = ctrl.state();
    near(s.motion.speed, 1 + t / (jr + jg) + jg / (jr + jg) * relative, 1e-10);
    near(s.generator_speed / 2, 1 + t / (jr + jg) - jr / (jr + jg) * relative, 1e-10);
    near(s.shaft_twist, (1 - std::cos(wn * t)) / (jr * wn * wn), 1e-10);
    auto copy = ctrl;
    ctrl.advance(.1, 2, .2);
    copy.advance(.1, 2, .2);
    near(ctrl.state().motion.azimuth, copy.state().motion.azimuth, 0);
    rejects([&] { copy.advance(.01, -1e20, 0); });
    check(copy.failed(), "Failed controller must not advance again");
    rejects([&] { copy.advance(.01, 1, 0); });
}
void kinematics(const Case &c) {
    BladeStructure blade(c);
    const ModalState modal{{.3, -.1, .2}, {0, 0, 0}};
    const auto station = blade.nodes[5];
    auto motion = [&](double t) {
        RotorKinematics k;
        k.azimuth = 1.1 * t + .5 * .2 * t * t;
        k.speed = 1.1 + .2 * t;
        k.acceleration = .2;
        k.pitch = .2 + .03 * t + .5 * .01 * t * t;
        k.pitch_rate = .03 + .01 * t;
        k.pitch_acceleration = .01;
        k.yaw = .1 + .02 * t + .5 * .005 * t * t;
        k.yaw_rate = .02 + .005 * t;
        k.yaw_acceleration = .005;
        blade.set_operation(k);
        return blade.motion(t, 0, modal, station);
    };
    const double t = .5, h = 1e-4;
    const auto minus = motion(t - h), center = motion(t), plus = motion(t + h);
    for (int i = 0; i < 3; ++i) {
        near((plus.position[i] - minus.position[i]) / (2 * h), center.velocity[i], 1e-6);
        near((plus.position[i] - 2 * center.position[i] + minus.position[i]) / (h * h),
             center.acceleration_bias[i], 1e-5);
    }
    Rotor rotor(c);
    const RotorState state{};
    const auto aerodynamic = rotor.evaluate(0, state);
    const auto loads = rotor.structural_loads(0, state, aerodynamic);
    double torque = 0;
    for (std::size_t b = 0; b < 3; ++b)
        for (std::size_t j = 0; j < loads[b].size(); ++j) {
            const auto motion = rotor.structure().motion(0, b, state[b], rotor.structure().nodes[j]);
            torque +=
                dot(cross(motion.position - rotor.structure().hub, loads[b][j].force) + loads[b][j].moment,
                    rotor.structure().shaft);
        }
    near(rotor.aerodynamic_torque(aerodynamic), torque, 1e-7);
}
void propagation() {
    using namespace aeroacoustics;
    Atmosphere air;
    check(air.absorption_db_per_m(1000) > 0, "Absorption must be positive");
    near(air.absorption_db_per_m(1000) * 1000, 4.98, .03); // rounded engineering table value, 20 C / 70% RH
    near(knife_edge_loss(0), 6.032852208563606, 1e-10);
    near(knife_edge_loss(-1), 0);
    Parameters p;
    p.freqlist = {1000};
    p.timod = 0;
    p.tbltemod = 1;
    Node n;
    n.aero_center = {0, 0, 20};
    n.section.speed = 40;
    // Proper orthonormal basis; edge sources differ only in y.
    n.global_to_local = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::vector<Node> nodes{n};
    std::vector<Vec3> observers{{100, 10, 0}};
    AcousticWorkspace free(p), ground(p), absorb(p);
    PropagationOptions o;
    o.absorption = false;
    o.ground = GroundModel::rigid;
    ground.set_propagation(
        std::make_shared<const OutdoorPropagation>(o, FrequencyBands::openfast_reference(p.freqlist)));
    auto reference = free.evaluate(nodes, observers), reflected = ground.evaluate(nodes, observers);
    for (std::size_t m = 0; m < mechanism_count; ++m)
        if (std::isfinite(reference[0][0][m][0]))
            near(reflected[0][0][m][0] - reference[0][0][m][0], 20 * std::log10(2.), 1e-10);
    o.ground = GroundModel::none;
    o.absorption = true;
    absorb.set_propagation(
        std::make_shared<const OutdoorPropagation>(o, FrequencyBands::openfast_reference(p.freqlist)));
    const auto attenuated = absorb.evaluate(nodes, observers);
    const auto paths =
        observe(observers[0], n.aero_center, n.global_to_local, n.section.chord, n.airfoil_reference);
    near(reference[0][0][index(Mechanism::trailing_pressure)][0] -
             attenuated[0][0][index(Mechanism::trailing_pressure)][0],
         air.absorption_db_per_m(1000) * paths.second.distance);
    o.absorption = false;
    o.screens = {{50, -100, 50, 100, 30}};
    absorb.set_propagation(
        std::make_shared<const OutdoorPropagation>(o, FrequencyBands::openfast_reference(p.freqlist)));
    const auto blocked = absorb.evaluate(nodes, observers);
    check(blocked[0][0][1][0] < reference[0][0][1][0] - 6, "Screen must attenuate blocked path");
    o.ground = GroundModel::rigid;
    absorb.set_propagation(
        std::make_shared<const OutdoorPropagation>(o, FrequencyBands::openfast_reference(p.freqlist)));
    near(absorb.evaluate(nodes, observers)[0][0][1][0] - blocked[0][0][1][0], 20 * std::log10(2.));
    o.normalized_impedance = {-1, 0};
    o.ground = GroundModel::impedance;
    rejects([&] { o.validate(); });
    rejects([&] {
        free.set_propagation(std::make_shared<const OutdoorPropagation>(
            PropagationOptions{}, FrequencyBands::openfast_reference({500})));
    });
}
void coupled(const Case &c, ControlConfig config) {
    struct SteadyAdapter final : WindField {
        SteadyWind field;
        explicit SteadyAdapter(SteadyWind w) : field(w) {}
        Vec3 at(double, const Vec3 &p) const override { return field.at(p); }
        const char *description() const noexcept override { return "steady adapter"; }
    };
    RunOptions baseline;
    baseline.duration = .1;
    Simulation original(c, baseline);
    baseline.solver.wind = std::make_shared<const SteadyAdapter>(c.wind);
    baseline.propagation = aeroacoustics::PropagationOptions{};
    baseline.propagation->absorption = false;
    Simulation degenerate(c, baseline);
    while (auto ref = original.next()) {
        auto test = degenerate.next();
        check(bool(test), "Missing degenerate frame");
        check(ref->generalized->q == test->generalized->q, "Steady field injection changed structure");
        if (ref->acoustics)
            check(ref->acoustics->power == test->acoustics->power, "No-loss propagation changed acoustics");
    }
    RunOptions o;
    o.duration = .2;
    o.solver.controller = config;
    o.solver.wind = std::make_shared<const GridWind>(std::array<std::vector<double>, 4>{{{0}, {0}, {0}, {0}}},
                                                     std::vector<Vec3>{{9, 1, 0}});
    o.propagation = aeroacoustics::PropagationOptions{};
    Simulation sim(c, o);
    sim.next();
    sim.next();
    auto checkpoint = sim.checkpoint();
    auto a = sim.next();
    auto state = *a->operation;
    auto q = a->generalized->q;
    sim.restore(checkpoint);
    auto b = sim.next();
    near(state.motion.azimuth, b->operation->motion.azimuth, 0);
    check(q == b->generalized->q, "Engineering restore changed structure");
    while (sim.next()) {
    };
    sim.summary();
    sim.reset();
    auto first = sim.next();
    near(first->operation->motion.azimuth, 0);
    const TurbineModel model(c);
    auto execute = [&] {
        Simulation run(model, o);
        double result = 0;
        while (auto frame = run.next())
            result += frame->operation->motion.azimuth;
        return result;
    };
    const double expected = execute();
    auto left = std::async(std::launch::async, execute), right = std::async(std::launch::async, execute);
    near(left.get(), expected, 0);
    near(right.get(), expected, 0);
}
} // namespace
int main(int argc, char **argv) {
    try {
        check(argc == 3, "Usage: engineering_probe CASE.fst controller.dat");
        Case c(argv[1]);
        auto config = ControlConfig::read(argv[2]);
        wind();
        drivetrain(config);
        kinematics(c);
        propagation();
        coupled(c, config);
        std::cout << "E1-E2 engineering analytic, propagation and checkpoint checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
