#include "turbine/tower.hpp"
#include "turbine/math.hpp"
#include <algorithm>
namespace turbine {
void TowerInfluence::validate() const {
    if (provenance.empty() || !std::isfinite(x) || !std::isfinite(y) || stations.size() < 2)
        throw std::invalid_argument("Invalid tower metadata/geometry");
    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto &s = stations[i];
        if (!std::isfinite(s.z) || !std::isfinite(s.diameter) || !std::isfinite(s.cd) || s.diameter <= 0 ||
            s.cd < 0 || s.cd > 3 || (i && s.z <= stations[i - 1].z))
            throw std::invalid_argument("Invalid tower station");
    }
}
std::shared_ptr<const TowerInfluence> TowerInfluence::read(const std::filesystem::path &p) {
    InputFile f(p);
    auto t = std::make_shared<TowerInfluence>();
    t->provenance = f.value("Provenance");
    t->x = f.number("X");
    t->y = f.number("Y");
    t->potential = f.flag("Potential");
    t->shadow = f.flag("Shadow");
    int n = f.integer("NumStations");
    if (n < 2 || n > 10000)
        throw std::invalid_argument("Tower station count");
    for (const auto &r : f.table_after("Stations", n, 3))
        t->stations.push_back({r[0], r[1], r[2]});
    t->validate();
    return t;
}
Vec3 TowerInfluence::apply(Vec3 wind, const Vec3 &p) const {
    for (double v : wind)
        if (!std::isfinite(v))
            throw std::invalid_argument("Invalid tower wind");
    for (double v : p)
        if (!std::isfinite(v))
            throw std::invalid_argument("Invalid tower position");
    if (!potential && !shadow)
        return wind;
    if (stations.size() < 2)
        throw std::invalid_argument("Uninitialized tower");
    const double z = std::clamp(p[2], stations.front().z, stations.back().z);
    auto upper = std::upper_bound(stations.begin(), stations.end(), z,
                                  [](double a, const TowerStation &b) { return a < b.z; });
    const auto j = std::clamp<std::size_t>(upper - stations.begin(), 1, stations.size() - 1);
    const auto &a = stations[j - 1];
    const auto &b = stations[j];
    const double q = (z - a.z) / (b.z - a.z);
    const double radius = .5 * (a.diameter + q * (b.diameter - a.diameter)), cd = a.cd + q * (b.cd - a.cd);
    const double height = std::abs(p[2] - z) / radius;
    if (height >= 1)
        return wind;
    const double speed = std::hypot(wind[0], wind[1]);
    if (speed < 1e-12)
        return wind;
    const Vec3 e{wind[0] / speed, wind[1] / speed, 0}, n{-e[1], e[0], 0}, offset{p[0] - x, p[1] - y, 0};
    const double clearance = std::hypot(norm(offset), p[2] - z) - radius;
    if (clearance <= 0)
        throw std::invalid_argument("Blade point intersects tower exclusion capsule");
    if (clearance > 40 * radius || clearance <= .02 * radius)
        return wind;
    const double taper = std::cos(.5 * pi * height);
    const double xx = dot(offset, e) / (radius * taper), yy = dot(offset, n) / (radius * taper),
                 rr = xx * xx + yy * yy;
    double du = potential ? (yy * yy - xx * xx) / (rr * rr) : 0,
           dv = potential ? -2 * xx * yy / (rr * rr) : 0;
    const double width = std::sqrt(std::sqrt(rr));
    if (shadow && xx > 0 && std::abs(yy) < width)
        du -= std::min(.5, cd / width * std::pow(std::cos(.5 * pi * yy / width), 2));
    return wind + speed * (du * e + dv * n);
}
} // namespace turbine
