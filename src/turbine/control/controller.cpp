#include "turbine/control.hpp"
namespace turbine {
ControlConfig ControlConfig::read(const std::filesystem::path &path) {
    InputFile f(path);
    ControlConfig c{};
    c.rotor_inertia = f.number("RotorInertia");
    c.generator_inertia = f.number("GeneratorInertia");
    c.gear_ratio = f.number("GearRatio");
    c.shaft_stiffness = f.number("ShaftStiffness");
    c.shaft_damping = f.number("ShaftDamping");
    c.optimal_torque_gain = f.number("OptimalTorqueGain");
    c.rated_power = f.number("RatedPower");
    c.efficiency = f.number("Efficiency");
    c.rated_rotor_speed = f.number("RatedRotorRPM") * 2 * pi / 60;
    c.torque_limit = f.number("TorqueLimit");
    c.torque_rate = f.number("TorqueRate");
    c.torque_time_constant = f.number("TorqueTimeConstant");
    c.pitch_min = f.number("PitchMinDeg") * deg;
    c.pitch_max = f.number("PitchMaxDeg") * deg;
    c.pitch_kp = f.number("PitchKp");
    c.pitch_ki = f.number("PitchKi");
    c.pitch_frequency = f.number("PitchFrequency");
    c.pitch_damping = f.number("PitchDamping");
    c.pitch_rate = f.number("PitchRateDeg") * deg;
    c.yaw_frequency = f.number("YawFrequency");
    c.yaw_damping = f.number("YawDamping");
    c.yaw_rate = f.number("YawRateDeg") * deg;
    c.yaw_deadband = f.number("YawDeadbandDeg") * deg;
    c.speed_filter_time = f.number("SpeedFilterTime");
    c.max_step = f.number("MaxStep");
    c.noise_start = f.number("NoiseStart");
    c.noise_speed_ratio = f.number("NoiseSpeedRatio");
    c.noise_power_ratio = f.number("NoisePowerRatio");
    c.noise_pitch = f.number("NoisePitchDeg") * deg;
    c.validate();
    return c;
}
void ControlConfig::validate() const {
    for (double v :
         {rotor_inertia,        generator_inertia, gear_ratio,        shaft_stiffness, optimal_torque_gain,
          rated_power,          efficiency,        rated_rotor_speed, torque_limit,    torque_rate,
          torque_time_constant, pitch_frequency,   pitch_damping,     pitch_rate,      yaw_frequency,
          yaw_damping,          yaw_rate,          speed_filter_time, max_step,        noise_speed_ratio,
          noise_power_ratio})
        if (!std::isfinite(v) || v <= 0)
            throw std::invalid_argument("Controller positive parameter invalid");
    for (double v : {shaft_damping, pitch_kp, pitch_ki, yaw_deadband, noise_start})
        if (!std::isfinite(v) || v < 0)
            throw std::invalid_argument("Controller nonnegative parameter invalid");
    if (!std::isfinite(pitch_min) || !std::isfinite(pitch_max) || !std::isfinite(noise_pitch) ||
        pitch_min < 0 || pitch_max > pi / 2 || pitch_max <= pitch_min || noise_pitch < pitch_min ||
        noise_pitch > pitch_max || efficiency > 1 || noise_speed_ratio > 1 || noise_power_ratio > 1)
        throw std::invalid_argument("Invalid controller limits or noise mode");
}
OperatingController::OperatingController(ControlConfig c, double speed, double pitch, double yaw)
    : config_(c) {
    c.validate();
    if (!std::isfinite(speed) || speed <= 0 || !std::isfinite(pitch) || pitch < c.pitch_min ||
        pitch > c.pitch_max || !std::isfinite(yaw))
        throw std::invalid_argument("Invalid controller initial state");
    x_[1] = speed;
    x_[2] = speed * c.gear_ratio;
    x_[5] = pitch;
    x_[9] = yaw;
    x_[8] = speed;
    x_[4] = std::min(c.torque_limit, c.optimal_torque_gain * x_[2] * x_[2]);
    x_[3] = c.gear_ratio * x_[4] / c.shaft_stiffness;
    x_[7] = pitch - c.pitch_min;
    publish(c.gear_ratio * x_[4], yaw);
}
OperatingController::State OperatingController::derivative(double t, const State &x, double aero,
                                                           double direction) const {
    const auto &c = config_;
    const bool noise = t >= c.noise_start;
    const double reference = c.rated_rotor_speed * (noise ? c.noise_speed_ratio : 1);
    const double floor = noise ? c.noise_pitch : c.pitch_min;
    const double error = x[8] - reference;
    const double demand = floor + c.pitch_kp * error + x[7];
    const double pitch_command = std::clamp(demand, floor, c.pitch_max);
    const double rated = c.rated_power * (noise ? c.noise_power_ratio : 1);
    const double torque_command =
        std::clamp(std::min(c.optimal_torque_gain * std::pow(std::max(0., x[8]) * c.gear_ratio, 2),
                            rated / (c.efficiency * std::max(x[2], .1))),
                   0., c.torque_limit);
    const double shaft = c.shaft_stiffness * x[3] + c.shaft_damping * (x[1] - x[2] / c.gear_ratio);
    State d{};
    d[0] = x[1];
    d[1] = (aero - shaft) / c.rotor_inertia;
    d[2] = (shaft / c.gear_ratio - x[4]) / c.generator_inertia;
    d[3] = x[1] - x[2] / c.gear_ratio;
    d[4] = std::clamp((torque_command - x[4]) / c.torque_time_constant, -c.torque_rate, c.torque_rate);
    d[5] = std::clamp(x[6], -c.pitch_rate, c.pitch_rate);
    d[6] = c.pitch_frequency * c.pitch_frequency * (pitch_command - x[5]) -
           2 * c.pitch_damping * c.pitch_frequency * x[6];
    if ((x[6] >= c.pitch_rate && d[6] > 0) || (x[6] <= -c.pitch_rate && d[6] < 0))
        d[6] = 0;
    d[7] = ((demand >= c.pitch_max && error > 0) || (demand <= floor && error < 0)) ? 0 : c.pitch_ki * error;
    d[8] = (x[2] / c.gear_ratio - x[8]) / c.speed_filter_time;
    const double yaw_error = std::remainder(direction - x[9], 2 * pi);
    d[9] = std::clamp(x[10], -c.yaw_rate, c.yaw_rate);
    d[10] = c.yaw_frequency * c.yaw_frequency * (std::abs(yaw_error) > c.yaw_deadband ? yaw_error : 0) -
            2 * c.yaw_damping * c.yaw_frequency * x[10];
    if ((x[10] >= c.yaw_rate && d[10] > 0) || (x[10] <= -c.yaw_rate && d[10] < 0))
        d[10] = 0;
    return d;
}
void OperatingController::publish(double torque, double direction) {
    const auto d = derivative(time_, x_, torque, direction);
    result_.motion = {x_[0], x_[1], d[1], x_[5], x_[6], d[6], x_[9], x_[10], d[10]};
    result_.generator_speed = x_[2];
    result_.shaft_twist = x_[3];
    result_.generator_torque = x_[4];
    result_.electrical_power = config_.efficiency * x_[4] * x_[2];
    result_.aerodynamic_torque = torque;
    result_.shaft_torque =
        config_.shaft_stiffness * x_[3] + config_.shaft_damping * (x_[1] - x_[2] / config_.gear_ratio);
    result_.noise_mode = time_ >= config_.noise_start;
    result_.speed_reference =
        config_.rated_rotor_speed * (result_.noise_mode ? config_.noise_speed_ratio : 1);
}
void OperatingController::advance(double dt, double torque, double direction) {
    if (failed_)
        throw std::logic_error("Controller failed; restore a valid copy or reconstruct before advancing");
    if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(torque) || !std::isfinite(direction))
        throw std::invalid_argument("Invalid controller input");
    const auto &c = config_;
    const double shaft_frequency = std::sqrt(
        c.shaft_stiffness * (1 / c.rotor_inertia + 1 / (c.generator_inertia * c.gear_ratio * c.gear_ratio)));
    const double shaft_decay =
        c.shaft_damping * (1 / c.rotor_inertia + 1 / (c.generator_inertia * c.gear_ratio * c.gear_ratio));
    const double step =
        std::min(c.max_step, .05 / std::max({shaft_frequency, shaft_decay, c.pitch_frequency,
                                             2 * c.pitch_damping * c.pitch_frequency, c.yaw_frequency,
                                             2 * c.yaw_damping * c.yaw_frequency, 1 / c.torque_time_constant,
                                             1 / c.speed_filter_time}));
    if (dt / step > 100000)
        throw std::runtime_error("Controller requires excessive substeps");
    const auto n = std::size_t(std::ceil(dt / step));
    const double h = dt / n;
    failed_ = true;
    auto trial = [](const State &x, const State &d, double h) {
        State y;
        for (std::size_t j = 0; j < x.size(); ++j)
            y[j] = x[j] + h * d[j];
        return y;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double t = time_ + i * h;
        auto a = derivative(t, x_, torque, direction),
             b = derivative(t + h / 2, trial(x_, a, h / 2), torque, direction);
        auto cc = derivative(t + h / 2, trial(x_, b, h / 2), torque, direction),
             d = derivative(t + h, trial(x_, cc, h), torque, direction);
        for (std::size_t j = 0; j < x_.size(); ++j) {
            x_[j] += h / 6 * (a[j] + 2 * b[j] + 2 * cc[j] + d[j]);
            if (!std::isfinite(x_[j]))
                throw std::runtime_error("Non-finite controller state");
        }
        x_[4] = std::clamp(x_[4], 0., c.torque_limit);
        x_[6] = std::clamp(x_[6], -c.pitch_rate, c.pitch_rate);
        x_[10] = std::clamp(x_[10], -c.yaw_rate, c.yaw_rate);
        if (x_[5] < c.pitch_min) {
            x_[5] = c.pitch_min;
            x_[6] = std::max(0., x_[6]);
        }
        if (x_[5] > c.pitch_max) {
            x_[5] = c.pitch_max;
            x_[6] = std::min(0., x_[6]);
        }
        if (x_[1] <= 0 || x_[2] <= 0)
            throw std::runtime_error("Rotor/generator stopped or reversed; startup/shutdown model required");
    }
    time_ += dt;
    publish(torque, direction);
    failed_ = false;
}
} // namespace turbine
