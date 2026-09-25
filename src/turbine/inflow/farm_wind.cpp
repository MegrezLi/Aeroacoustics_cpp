#include "turbine/farm.hpp"
#include <algorithm>
namespace turbine {
void WakeHistory::validate() const {
    if (!std::isfinite(radius) || radius <= 0 || samples.size() < 2)
        throw std::invalid_argument("Invalid wake history geometry/length");
    for (double x : hub)
        if (!std::isfinite(x))
            throw std::invalid_argument("Invalid wake hub");
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto &s = samples[i];
        if (!std::isfinite(s.time) || !std::isfinite(s.inflow_speed) || !std::isfinite(s.ct) ||
            s.inflow_speed <= 0 || s.ct < 0 || s.ct >= 1 || (i && s.time <= samples[i - 1].time))
            throw std::invalid_argument("Invalid wake history sample");
    }
}
WakeSample WakeHistory::at(double t) const {
    if (!std::isfinite(t) || samples.size() < 2 || t < samples.front().time - 1e-10 ||
        t > samples.back().time + 1e-10)
        throw std::out_of_range("Wake history extrapolation");
    t = std::clamp(t, samples.front().time, samples.back().time);
    auto it = std::upper_bound(samples.begin(), samples.end(), t,
                               [](double x, const WakeSample &s) { return x < s.time; });
    const auto j = std::clamp<std::size_t>(it - samples.begin(), 1, samples.size() - 1);
    const auto &a = samples[j - 1];
    const auto &b = samples[j];
    const double q = (t - a.time) / (b.time - a.time);
    return {t, a.inflow_speed + q * (b.inflow_speed - a.inflow_speed), a.ct + q * (b.ct - a.ct)};
}
double jensen_deficit(double ct, double u, double r, double k, double x, double transverse) {
    for (double v : {ct, u, r, k, x, transverse})
        if (!std::isfinite(v))
            throw std::invalid_argument("Non-finite Jensen input");
    if (ct < 0 || ct >= 1 || u <= 0 || r <= 0 || k <= 0 || transverse < 0)
        throw std::invalid_argument("Invalid Jensen input");
    if (x <= 0 || transverse > r + k * x)
        return 0;
    return u * (ct / (1 + std::sqrt(1 - ct))) * std::pow(r / (r + k * x), 2);
}
FarmWind::FarmWind(SteadyWind a, Vec3 origin, double k,
                   std::vector<std::shared_ptr<const WakeHistory>> upstream)
    : ambient_(a), origin_(origin), direction_{std::cos(a.propagation), -std::sin(a.propagation), 0},
      expansion_(k), upstream_(std::move(upstream)) {
    if (!std::isfinite(k) || k <= 0 || !std::isfinite(a.speed) || a.speed <= 0 || a.upflow != 0)
        throw std::invalid_argument("Farm requires positive steady horizontal wind/expansion");
    for (double v : origin_)
        if (!std::isfinite(v))
            throw std::invalid_argument("Invalid farm origin");
    for (const auto &p : upstream_) {
        if (!p)
            throw std::invalid_argument("Null wake history");
        p->validate();
    }
}
Vec3 FarmWind::at(double t, const Vec3 &local) const {
    if (!std::isfinite(t))
        throw std::invalid_argument("Invalid farm query time");
    for (double v : local)
        if (!std::isfinite(v))
            throw std::invalid_argument("Invalid farm query position");
    const Vec3 p = local + origin_;
    const Vec3 base = ambient_.at(p);
    double deficit2 = 0;
    for (const auto &w : upstream_) {
        const auto state = w->at(t);
        const Vec3 d = local + (origin_ - w->hub);
        const double x = dot(d, direction_);
        const double lateral = norm(d - x * direction_);
        const double du = jensen_deficit(state.ct, state.inflow_speed, w->radius, expansion_, x, lateral);
        deficit2 += du * du;
    }
    const double deficit = std::sqrt(deficit2);
    if (!std::isfinite(deficit) || deficit >= dot(base, direction_))
        throw std::invalid_argument("Overlapping Jensen wakes exhaust local wind; case outside model domain");
    return base - deficit * direction_;
}
} // namespace turbine
