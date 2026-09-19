#pragma once
#include "turbine/math.hpp"
namespace turbine {
struct RotorKinematics {
    double azimuth = 0, speed = 0, acceleration = 0;
    double pitch = 0, pitch_rate = 0, pitch_acceleration = 0;
    double yaw = 0, yaw_rate = 0, yaw_acceleration = 0;
};
struct ControlConfig {
    double rotor_inertia, generator_inertia, gear_ratio, shaft_stiffness, shaft_damping;
    double optimal_torque_gain, rated_power, efficiency, rated_rotor_speed;
    double torque_limit, torque_rate, torque_time_constant;
    double pitch_min, pitch_max, pitch_kp, pitch_ki, pitch_frequency, pitch_damping, pitch_rate;
    double yaw_frequency, yaw_damping, yaw_rate, yaw_deadband;
    double speed_filter_time, max_step;
    double noise_start, noise_speed_ratio, noise_power_ratio, noise_pitch;
    static ControlConfig read(const std::filesystem::path &);
    void validate() const;
};
struct ControlState {
    RotorKinematics motion;
    double generator_speed = 0, shaft_twist = 0, generator_torque = 0, electrical_power = 0;
    double aerodynamic_torque = 0, shaft_torque = 0, speed_reference = 0;
    bool noise_mode = false;
};
// Two-mass drivetrain with a first-order generator torque response, collective
// PI pitch control, and second-order pitch/yaw servos. Aerodynamic torque is held
// over a structural step; internal RK4 steps resolve shaft/actuator dynamics.
class OperatingController {
    ControlConfig config_;
    // azimuth, rotor speed, generator speed, shaft twist, torque, pitch, pitch rate,
    // pitch integral, filtered rotor-equivalent speed, yaw, yaw rate
    using State = std::array<double, 11>;
    State x_{};
    double time_ = 0;
    bool failed_ = false;
    ControlState result_;
    State derivative(double, const State &, double, double) const;
    void publish(double torque, double direction);

  public:
    OperatingController(ControlConfig, double rotor_speed, double pitch, double yaw);
    const ControlState &state() const noexcept { return result_; }
    const ControlConfig &config() const noexcept { return config_; }
    bool failed() const noexcept { return failed_; }
    void advance(double dt, double aerodynamic_torque, double wind_direction);
};
} // namespace turbine
