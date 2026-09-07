#pragma once
#include "turbine/input.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace turbine {
using Matrix3 = std::array<Vec3, 3>;
inline Vec3 operator+(Vec3 a, const Vec3 &b) {
    for (int i = 0; i < 3; ++i)
        a[i] += b[i];
    return a;
}
inline Vec3 operator-(Vec3 a, const Vec3 &b) {
    for (int i = 0; i < 3; ++i)
        a[i] -= b[i];
    return a;
}
inline Vec3 operator-(Vec3 a) {
    for (auto &x : a)
        x = -x;
    return a;
}
inline Vec3 operator*(double s, Vec3 a) {
    for (auto &x : a)
        x *= s;
    return a;
}
inline Vec3 operator*(Vec3 a, double s) { return s * a; }
inline Vec3 operator/(Vec3 a, double s) { return (1 / s) * a; }
inline double dot(const Vec3 &a, const Vec3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const Vec3 &a) { return std::sqrt(dot(a, a)); }
inline Vec3 unit(const Vec3 &a) {
    double n = norm(a);
    if (n < 1e-14)
        throw std::runtime_error("Zero-length direction");
    return a / n;
}
inline Vec3 multiply(const Matrix3 &m, const Vec3 &v) { return {dot(m[0], v), dot(m[1], v), dot(m[2], v)}; }
inline Matrix3 transpose(const Matrix3 &a) {
    Matrix3 b{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            b[i][j] = a[j][i];
    return b;
}
inline Matrix3 multiply(const Matrix3 &a, const Matrix3 &b) {
    Matrix3 c{};
    auto bt = transpose(b);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c[i][j] = dot(a[i], bt[j]);
    return c;
}
inline Matrix3 euler_matrix(const Vec3 &a) {
    const double cx = std::cos(a[0]), sx = std::sin(a[0]), cy = std::cos(a[1]), sy = std::sin(a[1]),
                 cz = std::cos(a[2]), sz = std::sin(a[2]);
    return {{{cy * cz, cx * sz + sx * sy * cz, sx * sz - cx * sy * cz},
             {-cy * sz, cx * cz - sx * sy * sz, sx * cz + cx * sy * sz},
             {sy, -sx * cy, cx * cy}}};
}
inline Vec3 euler_angles(const Matrix3 &m) {
    double cy = std::hypot(m[0][0], m[1][0]);
    if (cy < 1e-14)
        return {std::atan2(m[1][2], m[1][1]), std::atan2(m[2][0], cy), 0};
    const double z = std::atan2(-m[1][0], m[0][0]), cz = std::cos(z), sz = std::sin(z);
    cy = std::copysign(cy, std::abs(cz) < 1e-14 ? -m[1][0] / sz : m[0][0] / cz);
    return {std::atan2(sz * m[0][2] + cz * m[1][2], sz * m[0][1] + cz * m[1][1]), std::atan2(m[2][0], cy), z};
}
inline Vec3 rotation_log(const Matrix3 &m) {
    const double c = std::clamp((m[0][0] + m[1][1] + m[2][2] - 1) / 2, -1.0, 1.0), t = std::acos(c);
    const Vec3 skew{m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]};
    if (t < 1e-14)
        return {};
    if (t <= 3.1)
        return (t / (2 * std::sin(t))) * skew;
    Vec3 v{1 + m[0][0] - m[1][1] - m[2][2], 1 - m[0][0] + m[1][1] - m[2][2], 1 - m[0][0] - m[1][1] + m[2][2]};
    int k = 0;
    for (int i = 1; i < 3; ++i)
        if (std::abs(v[i]) > std::abs(v[k]))
            k = i;
    const double d = std::sqrt(std::abs(v[k] * 2 * (1 - c))) / t;
    for (int i = 0; i < 3; ++i)
        if (i != k)
            v[i] = m[i][k] + m[k][i];
    v = v / d;
    if (std::abs(t - pi) > 1e-14) {
        k = 0;
        for (int i = 1; i < 3; ++i)
            if (std::abs(skew[i]) > std::abs(skew[k]))
                k = i;
        if (std::signbit(skew[k]) != std::signbit(v[k]))
            v = -v;
    }
    return v;
}
inline Matrix3 rotation_exp(const Vec3 &v) {
    const double t = norm(v), a = t < 1e-8 ? 1 - t * t / 6 : std::sin(t) / t,
                 b = t < 1e-8 ? .5 - t * t / 24 : (1 - std::cos(t)) / (t * t);
    const double x = v[0], y = v[1], z = v[2];
    return {{{1 - b * (y * y + z * z), b * x * y + a * z, b * x * z - a * y},
             {b * x * y - a * z, 1 - b * (x * x + z * z), b * y * z + a * x},
             {b * x * z + a * y, b * y * z - a * x, 1 - b * (x * x + y * y)}}};
}
inline Matrix3 interpolate_rotation(const Matrix3 &a, const Matrix3 &b, double w) {
    if (w < 1e-14)
        return a;
    if (w > 1 - 1e-14)
        return b;
    const Vec3 va = rotation_log(a);
    Vec3 vb = rotation_log(b);
    if (norm(vb) > 1e-14) {
        const Vec3 period = (2 * pi / norm(vb)) * vb;
        const double shift = std::round(dot(va - vb, period) / dot(period, period));
        vb = vb + shift * period;
    }
    return rotation_exp((1 - w) * va + w * vb);
}
inline Matrix3 small_rotation(const Vec3 &a) {
    const double w2 = 1 - dot(a, a) / 4;
    if (w2 < 0)
        throw std::runtime_error("Small-rotation domain exceeded");
    const double w = std::sqrt(w2);
    const double x = a[0], y = a[1], z = a[2];
    return {{{1 - (y * y + z * z) / 2, x * y / 2 + z * w, x * z / 2 - y * w},
             {x * y / 2 - z * w, 1 - (x * x + z * z) / 2, y * z / 2 + x * w},
             {x * z / 2 + y * w, y * z / 2 - x * w, 1 - (x * x + y * y) / 2}}};
}
inline Vec3 solve3(Matrix3 a, Vec3 b) {
    for (int i = 0; i < 3; ++i) {
        int pivot = i;
        for (int j = i + 1; j < 3; ++j)
            if (std::abs(a[j][i]) > std::abs(a[pivot][i]))
                pivot = j;
        if (std::abs(a[pivot][i]) < 1e-16)
            throw std::runtime_error("Singular structural mass matrix");
        std::swap(a[i], a[pivot]);
        std::swap(b[i], b[pivot]);
        for (int j = i + 1; j < 3; ++j) {
            const double f = a[j][i] / a[i][i];
            for (int k = i; k < 3; ++k)
                a[j][k] -= f * a[i][k];
            b[j] -= f * b[i];
        }
    }
    Vec3 x{};
    for (int i = 2; i >= 0; --i) {
        double r = b[i];
        for (int j = i + 1; j < 3; ++j)
            r -= a[i][j] * x[j];
        x[i] = r / a[i][i];
    }
    return x;
}
} // namespace turbine
