#include "turbine/structure.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace turbine {
namespace {
double shape(const std::array<double, 5> &c, double x, double L, int derivative) {
    double y = 0;
    for (int k = 2; k <= 6; ++k) {
        double factor = 1;
        for (int j = 0; j < derivative; ++j)
            factor *= k - j;
        y += c[k - 2] * factor * std::pow(x, k - derivative) / std::pow(L, derivative);
    }
    return y;
}
} // namespace
BladeStructure::BladeStructure(const Case &c)
    : length(c.structure.number("TipRad") - c.structure.number("HubRad")),
      hub_radius(c.structure.number("HubRad")), omega(c.structure.number("RotSpeed") * 2 * pi / 60),
      initial_azimuth((c.structure.number("Azimuth") - c.structure.number("AzimB1Up")) * deg - pi / 2),
      gravity_(c.gravity), tilt_(c.structure.number("ShftTilt") * deg),
      yaw_(c.structure.number("NacYaw") * deg) {
    if (c.structure.number("OoPDefl") != 0 || c.structure.number("IPDefl") != 0)
        throw std::runtime_error("Initial blade deflection initialization is not implemented yet");
    shaft = {std::cos(tilt_) * std::cos(yaw_), std::cos(tilt_) * std::sin(yaw_), std::sin(tilt_)};
    hub = c.structure.number("OverHang") * shaft +
          Vec3{0, 0, c.structure.number("TowerHt") + c.structure.number("Twr2Shft")};
    overhang_ = c.structure.number("OverHang");
    yaw_pivot_ = {0, 0, c.structure.number("TowerHt") + c.structure.number("Twr2Shft")};
    for (int b = 0; b < 3; ++b) {
        const auto suffix = "(" + std::to_string(b + 1) + ")";
        cone[b] = c.structure.number("PreCone" + suffix) * deg;
        pitch[b] = c.structure.number("BlPitch" + suffix) * deg;
        tip_mass[b] = c.structure.number("TipMass" + suffix);
    }
    const auto &f = c.blade_structure;
    const int count = c.structure.integer("BldNodes");
    if (count < 2 || length <= 0)
        throw std::runtime_error("Invalid blade structural discretization");
    std::array<std::array<double, 5>, 3> polynomial{};
    const std::array<std::string, 3> keys{{"BldFl1Sh", "BldFl2Sh", "BldEdgSh"}};
    for (int m = 0; m < 3; ++m)
        for (int k = 2; k <= 6; ++k)
            polynomial[m][k - 2] = f.number(keys[m] + "(" + std::to_string(k) + ")");
    std::vector<std::array<double, 5>> table;
    bool found = false;
    for (const auto &row : f.rows) {
        if (!row.empty() && row[0] == "BlFract") {
            found = true;
            continue;
        }
        if (!found || row.size() < 5)
            continue;
        try {
            std::array<double, 5> v{};
            for (int i = 0; i < 5; ++i)
                v[i] = std::stod(row[i]);
            table.push_back(v);
            if (table.size() == static_cast<std::size_t>(f.integer("NBlInpSt")))
                break;
        } catch (const std::invalid_argument &) {
        }
    }
    if (table.size() != static_cast<std::size_t>(f.integer("NBlInpSt")))
        throw std::runtime_error("Incomplete structural property table");
    auto property = [&](double s, int col) {
        if (s <= table.front()[0])
            return table.front()[col];
        if (s >= table.back()[0])
            return table.back()[col];
        auto it =
            std::upper_bound(table.begin(), table.end(), s, [](double x, const auto &r) { return x < r[0]; });
        auto l = it - 1;
        double w = (s - (*l)[0]) / ((*it)[0] - (*l)[0]);
        return (*l)[col] + w * ((*it)[col] - (*l)[col]);
    };
    nodes.resize(count + 2);
    nodes.front().twist = table.front()[1] * deg;
    nodes.back().twist = table.back()[1] * deg;
    nodes.back().span = length;
    const double dr = length / count;
    Vec3 old_flap_half{}, old_edge_half{}, old_flap_slope_half{}, old_edge_slope_half{};
    Matrix3 old_axial_half{};
    for (int j = 1; j <= count; ++j) {
        auto &s = nodes[j];
        const auto &prev = nodes[j - 1];
        s.span = (j - .5) * dr;
        s.width = dr;
        double fraction = s.span / length;
        s.mass = property(fraction, 2) * f.number("AdjBlMs") * dr;
        s.twist = property(fraction, 1) * deg;
        Vec3 curvature{};
        for (int m = 0; m < 3; ++m) {
            const double v = shape(polynomial[m], fraction, length, 0);
            modal_mass[m] += s.mass * v * v;
            curvature[m] = shape(polynomial[m], fraction, length, 2);
        }
        const double flapst = property(fraction, 3) * f.number("AdjFlSt") * dr,
                     edgest = property(fraction, 4) * f.number("AdjEdSt") * dr;
        for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
                stiffness[m][n] += flapst * curvature[m] * curvature[n];
        stiffness[2][2] += edgest * curvature[2] * curvature[2];
        const double co = std::cos(s.twist), si = std::sin(s.twist);
        const Vec3 cflap{curvature[0] * co, curvature[1] * co, curvature[2] * si},
            cedge{-curvature[0] * si, -curvature[1] * si, curvature[2] * co};
        const Vec3 flap_slope_half = .5 * dr * cflap, edge_slope_half = .5 * dr * cedge;
        s.flap_slope = flap_slope_half + prev.flap_slope + old_flap_slope_half;
        s.edge_slope = edge_slope_half + prev.edge_slope + old_edge_slope_half;
        const Vec3 flap_half = .5 * dr * s.flap_slope, edge_half = .5 * dr * s.edge_slope;
        s.flap = flap_half + prev.flap + old_flap_half;
        s.edge = edge_half + prev.edge + old_edge_half;
        for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n) {
                const double half =
                    .5 * dr * (s.flap_slope[m] * s.flap_slope[n] + s.edge_slope[m] * s.edge_slope[n]);
                s.axial[m][n] = half + prev.axial[m][n] + old_axial_half[m][n];
                old_axial_half[m][n] = half;
            }
        old_flap_slope_half = flap_slope_half;
        old_edge_slope_half = edge_slope_half;
        old_flap_half = flap_half;
        old_edge_half = edge_half;
    }
    auto &tip = nodes.back();
    const auto &last = nodes[count];
    tip.flap = last.flap + old_flap_half;
    tip.edge = last.edge + old_edge_half;
    tip.flap_slope = last.flap_slope + old_flap_slope_half;
    tip.edge_slope = last.edge_slope + old_edge_slope_half;
    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            tip.axial[m][n] = last.axial[m][n] + old_axial_half[m][n];
    const Vec3 tune{f.number("FlStTunr1"), f.number("FlStTunr2"), 1},
        damp{f.number("BldFlDmp1"), f.number("BldFlDmp2"), f.number("BldEdDmp1")};
    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            stiffness[m][n] *= std::sqrt(tune[m] * tune[n]);
    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            damping[m][n] = .02 * damp[n] * stiffness[m][n] / std::sqrt(stiffness[n][n] / modal_mass[n]);
}
void BladeStructure::set_operation(const RotorKinematics &op) {
    for (double v : {op.azimuth, op.speed, op.acceleration, op.pitch, op.pitch_rate, op.pitch_acceleration,
                     op.yaw, op.yaw_rate, op.yaw_acceleration})
        if (!std::isfinite(v))
            throw std::invalid_argument("Non-finite rotor kinematics");
    operation_ = op;
    omega = op.speed;
    pitch.fill(op.pitch);
    yaw_ = op.yaw;
    shaft = {std::cos(tilt_) * std::cos(yaw_), std::cos(tilt_) * std::sin(yaw_), std::sin(tilt_)};
    hub = overhang_ * shaft + yaw_pivot_;
}
Matrix3 BladeStructure::blade_basis(double time, std::size_t blade) const {
    const Vec3 c2{-std::sin(tilt_) * std::cos(yaw_), -std::sin(tilt_) * std::sin(yaw_), std::cos(tilt_)},
        c3{std::sin(yaw_), -std::cos(yaw_), 0};
    const double a = initial_azimuth + (operation_ ? operation_->azimuth : omega * time) + 2 * pi * blade / 3;
    const Vec3 g2 = std::cos(a) * c2 + std::sin(a) * c3, g3 = -std::sin(a) * c2 + std::cos(a) * c3;
    const Vec3 i1 = std::cos(cone[blade]) * shaft - std::sin(cone[blade]) * g3,
               i3 = std::sin(cone[blade]) * shaft + std::cos(cone[blade]) * g3;
    return {{std::cos(pitch[blade]) * i1 - std::sin(pitch[blade]) * g2,
             std::sin(pitch[blade]) * i1 + std::cos(pitch[blade]) * g2, i3}};
}
Motion BladeStructure::motion(double time, std::size_t blade, const ModalState &state,
                              const StructuralStation &s) const {
    return motion(blade_basis(time, blade), state, s);
}
void BladeStructure::motions_into(const Matrix3 &basis, const ModalState &state,
                                  std::vector<Motion> &result) const {
    result.resize(nodes.size());
    for (std::size_t j = 0; j < nodes.size(); ++j)
        result[j] = motion(basis, state, nodes[j]);
}
Motion BladeStructure::motion(const Matrix3 &b, const ModalState &state, const StructuralStation &s) const {
    const auto &q = state.q;
    const auto &qd = state.qd;
    Motion y;
    const Vec3 axialq = multiply(s.axial, q), axialqd = multiply(s.axial, qd);
    const Vec3 relative =
        dot(s.flap, q) * b[0] + dot(s.edge, q) * b[1] + (hub_radius + s.span - .5 * dot(q, axialq)) * b[2];
    y.position = hub + relative;
    Vec3 flexible_velocity{};
    for (int m = 0; m < 3; ++m) {
        y.partial_velocity[m] = s.flap[m] * b[0] + s.edge[m] * b[1] - axialq[m] * b[2];
        y.partial_angular[m] = -s.edge_slope[m] * b[0] + s.flap_slope[m] * b[1];
        flexible_velocity = flexible_velocity + qd[m] * y.partial_velocity[m];
    }
    Vec3 om = omega * shaft, angular_acceleration{}, hub_velocity{}, hub_acceleration{};
    if (operation_) {
        const auto &op = *operation_;
        const Vec3 z{0, 0, 1}, yaw_velocity = op.yaw_rate * z;
        const Vec3 carrier = yaw_velocity + omega * shaft;
        om = carrier - op.pitch_rate * b[2];
        angular_acceleration = op.yaw_acceleration * z + op.acceleration * shaft +
                               omega * cross(yaw_velocity, shaft) - op.pitch_acceleration * b[2] -
                               op.pitch_rate * cross(carrier, b[2]);
        hub_velocity = cross(yaw_velocity, hub - yaw_pivot_);
        hub_acceleration =
            cross(op.yaw_acceleration * z, hub - yaw_pivot_) + cross(yaw_velocity, hub_velocity);
    }
    y.velocity = cross(om, relative) + flexible_velocity;
    y.acceleration_bias =
        cross(om, cross(om, relative)) + 2 * cross(om, flexible_velocity) - dot(qd, axialqd) * b[2];
    if (operation_) {
        y.velocity = y.velocity + hub_velocity;
        y.acceleration_bias = y.acceleration_bias + hub_acceleration + cross(angular_acceleration, relative);
    }
    y.angular_velocity = om;
    for (int m = 0; m < 3; ++m)
        y.angular_velocity = y.angular_velocity + qd[m] * y.partial_angular[m];
    const double co = std::cos(s.twist), si = std::sin(s.twist), oop = dot(s.flap_slope, q),
                 ip = -dot(s.edge_slope, q);
    const Matrix3 local{{co * b[0] - si * b[1], si * b[0] + co * b[1], b[2]}};
    y.orientation = multiply(small_rotation({co * ip - si * oop, si * ip + co * oop, 0}), local);
    return y;
}
Vec3 BladeStructure::acceleration(double time, std::size_t blade, const ModalState &s,
                                  const std::vector<PointLoad> &loads) const {
    std::vector<Motion> motions;
    motions_into(blade_basis(time, blade), s, motions);
    return acceleration(blade, s, motions, loads);
}
Vec3 BladeStructure::acceleration(std::size_t blade, const ModalState &s, const std::vector<Motion> &motions,
                                  const std::vector<PointLoad> &loads) const {
    if (blade >= tip_mass.size() || motions.size() != nodes.size() || loads.size() != nodes.size())
        throw std::runtime_error("Structural load count mismatch");
    Matrix3 mass{};
    Vec3 force{};
    for (std::size_t j = 1; j < nodes.size(); ++j) {
        const auto &motion_j = motions[j];
        const double m = j + 1 == nodes.size() ? tip_mass[blade] : nodes[j].mass;
        const Vec3 f = loads[j].force - m * (Vec3{0, 0, gravity_} + motion_j.acceleration_bias);
        for (int i = 0; i < 3; ++i) {
            force[i] +=
                dot(motion_j.partial_velocity[i], f) + dot(motion_j.partial_angular[i], loads[j].moment);
            for (int k = 0; k < 3; ++k)
                mass[i][k] += m * dot(motion_j.partial_velocity[i], motion_j.partial_velocity[k]);
        }
    }
    force = force - multiply(stiffness, s.q) - multiply(damping, s.qd);
    return solve3(mass, force);
}
} // namespace turbine
