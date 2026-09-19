#include "acoustic_levels.hpp"
#include "propagation.hpp"
#include <algorithm>
namespace aeroacoustics {
namespace {
constexpr double pi = 3.14159265358979323846;
double distance(const Vec3 &a, const Vec3 &b) {
    return std::hypot(std::hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}
Vec3 source(const Node &n, bool leading) {
    Vec3 p = n.aero_center;
    for (int i = 0; i < 3; ++i)
        p[i] += n.section.chord * (-n.airfoil_reference[1] * n.global_to_local[i] +
                                   ((leading ? 0. : 1.) - n.airfoil_reference[0]) * n.global_to_local[3 + i]);
    return p;
}
} // namespace
void Atmosphere::validate() const {
    if (!std::isfinite(temperature_k) || temperature_k < 253.15 || temperature_k > 323.15 ||
        !std::isfinite(relative_humidity_percent) || relative_humidity_percent < 10 ||
        relative_humidity_percent > 100 || !std::isfinite(pressure_pa) || pressure_pa < 50000 ||
        pressure_pa > 120000)
        throw std::invalid_argument("Atmosphere requires -20..50 C, 10..100 percent RH, 50..120 kPa");
}
double Atmosphere::absorption_db_per_m(double f) const {
    validate();
    if (!std::isfinite(f) || f <= 0)
        throw std::invalid_argument("Invalid absorption frequency");
    const double tr = temperature_k / 293.15, pr = pressure_pa / 101325.;
    const double h = relative_humidity_percent *
                     std::pow(10., -6.8346 * std::pow(273.16 / temperature_k, 1.261) + 4.6151) / pr;
    const double fo = pr * (24 + 4.04e4 * h * (.02 + h) / (.391 + h));
    const double fn = pr * std::pow(tr, -.5) * (9 + 280 * h * std::exp(-4.170 * (std::pow(tr, -1. / 3) - 1)));
    return 8.686 * f * f *
           (1.84e-11 / pr * std::sqrt(tr) +
            std::pow(tr, -2.5) * (.01275 * std::exp(-2239.1 / temperature_k) / (fo + f * f / fo) +
                                  .1068 * std::exp(-3352. / temperature_k) / (fn + f * f / fn)));
}
void PropagationOptions::validate() const {
    atmosphere.validate();
    if (!std::isfinite(ground_z) || !std::isfinite(sound_speed) || sound_speed <= 0 ||
        (ground != GroundModel::none && ground != GroundModel::rigid && ground != GroundModel::impedance))
        throw std::invalid_argument("Invalid propagation geometry/model");
    if (ground == GroundModel::impedance &&
        (!std::isfinite(normalized_impedance.real()) || !std::isfinite(normalized_impedance.imag()) ||
         normalized_impedance.real() <= 0))
        throw std::invalid_argument("Ground impedance must be finite and passive (Re Z > 0)");
    for (const auto &s : screens) {
        for (double x : {s.x1, s.y1, s.x2, s.y2, s.top_z})
            if (!std::isfinite(x))
                throw std::invalid_argument("Non-finite screen");
        if (std::hypot(s.x2 - s.x1, s.y2 - s.y1) <= 0 || s.top_z <= ground_z)
            throw std::invalid_argument("Invalid screen dimensions");
    }
}
double knife_edge_loss(double v) {
    if (!std::isfinite(v))
        throw std::invalid_argument("Non-finite Fresnel parameter");
    return v <= -.78 ? 0 : 6.9 + 20 * std::log10(std::hypot(v - .1, 1.) + v - .1);
}
OutdoorPropagation::OutdoorPropagation(PropagationOptions o, FrequencyBands b)
    : options_(std::move(o)), bands_(std::move(b)) {
    options_.validate();
    for (const auto &band : bands_.values())
        absorption_.push_back(options_.absorption ? options_.atmosphere.absorption_db_per_m(band.center_hz)
                                                  : 0.);
}
Vec3 OutdoorPropagation::image_observer(Vec3 p) const {
    p[2] = 2 * options_.ground_z - p[2];
    return p;
}
double OutdoorPropagation::screen_loss(const Vec3 &s, const Vec3 &r, double f) const {
    const double dx = r[0] - s[0], dy = r[1] - s[1], horizontal = std::hypot(dx, dy);
    if (horizontal == 0)
        return 0;
    double loss = 0;
    for (const auto &edge : options_.screens) {
        const double ex = edge.x2 - edge.x1, ey = edge.y2 - edge.y1;
        const double det = dx * ey - dy * ex;
        if (std::abs(det) < 1e-12 * horizontal * std::hypot(ex, ey))
            continue;
        const double a = ((edge.x1 - s[0]) * ey - (edge.y1 - s[1]) * ex) / det;
        const double b = ((edge.x1 - s[0]) * dy - (edge.y1 - s[1]) * dx) / det;
        if (a <= 0 || a >= 1 || b < 0 || b > 1)
            continue;
        const double h = edge.top_z - (s[2] + a * (r[2] - s[2]));
        const double direct = distance(s, r);
        const double v =
            h * (horizontal / direct) * std::sqrt(2 * f / options_.sound_speed / (direct * a * (1 - a)));
        loss = std::max(loss, knife_edge_loss(v));
    }
    return loss;
}
void OutdoorPropagation::apply(const Node &node, const Vec3 &receiver, Mechanisms &direct,
                               const Mechanisms *image) const {
    for (double x : receiver)
        if (!std::isfinite(x))
            throw std::invalid_argument("Invalid propagation receiver");
    if (receiver[2] < options_.ground_z)
        throw std::invalid_argument("Receiver below ground plane");
    if (has_reflection() && !image)
        throw std::invalid_argument("Missing reflected source spectrum");
    for (std::size_t m = 0; m < mechanism_count; ++m) {
        if (direct[m].size() != bands_.size() || (image && (*image)[m].size() != bands_.size()))
            throw std::invalid_argument("Propagation frequency shape mismatch");
        const auto s = source(node, m == index(Mechanism::inflow));
        if (s[2] < options_.ground_z)
            throw std::invalid_argument("Source below ground plane");
        const double rd = distance(s, receiver), rr = distance(s, image_observer(receiver));
        if (rd <= 0 || rr <= 0)
            throw std::invalid_argument("Receiver coincides with source");
        for (std::size_t i = 0; i < bands_.size(); ++i) {
            const auto &band = bands_.values()[i];
            const double diffraction = screen_loss(s, receiver, band.center_hz);
            const double direct_db = direct[m][i] - absorption_[i] * rd - diffraction;
            if (!has_reflection()) {
                direct[m][i] = direct_db;
                continue;
            }
            std::complex<double> reflection{1, 0};
            if (options_.ground == GroundModel::impedance) {
                const double cosine = (s[2] + receiver[2] - 2 * options_.ground_z) / rr;
                const auto z = options_.normalized_impedance * cosine;
                reflection = (z - 1.) / (z + 1.);
            }
            const double height_sum = s[2] + receiver[2] - 2 * options_.ground_z;
            const double fraction = height_sum > 0 ? (s[2] - options_.ground_z) / height_sum : .5;
            const Vec3 reflection_point{s[0] + fraction * (receiver[0] - s[0]),
                                        s[1] + fraction * (receiver[1] - s[1]), options_.ground_z};
            const double reflection_loss = std::max(screen_loss(s, reflection_point, band.center_hz),
                                                    screen_loss(reflection_point, receiver, band.center_hz));
            const double image_db = (*image)[m][i] - absorption_[i] * rr - reflection_loss;
            const double maximum = std::max(direct_db, image_db);
            if (maximum == -INFINITY) {
                direct[m][i] = -INFINITY;
                continue;
            }
            // Stable linear sum. Flat within-band source density; integrate path
            // phase analytically, rather than aliasing long-path interference.
            const double d = relative_power(direct_db - maximum), r = relative_power(image_db - maximum);
            const double delay = (rr - rd) / options_.sound_speed;
            const double x = pi * (band.upper_hz - band.lower_hz) * delay;
            const double sinc = std::abs(x) < 1e-10 ? 1 : std::sin(x) / x;
            const auto phase = std::polar(1., -2 * pi * .5 * (band.lower_hz + band.upper_hz) * delay);
            const double energy =
                std::max(0., d + std::norm(reflection) * r +
                                 2 * std::sqrt(d * r) * std::real(reflection * phase) * sinc);
            direct[m][i] = energy == 0 ? -INFINITY : maximum + 10 * std::log10(energy);
        }
    }
}
} // namespace aeroacoustics
