#include "turbine/unsteady.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace turbine {
namespace {
double decay(double dt, double T, double old, double x, double oldx) {
    return old * std::exp(-dt / T) + (x - oldx) * std::exp(-.5 * dt / T);
}
double sq(double x) { return x * x; }
double blend(double x, double lo, double hi) {
    if (x <= lo)
        return 0;
    if (x >= hi)
        return 1;
    return .5 * (1 - std::cos(pi * (x - lo) / (hi - lo)));
}
bool small_single(double x) { return std::abs(x) < 1e-6; }
} // namespace
UnsteadyAirfoil::UnsteadyAirfoil(const Airfoil &af, double c, double dt, double sound)
    : af_(&af), chord_(c), dt_(dt), sound_(sound) {
    const auto &f = af.input;
    p_ = {f.number("alpha0", 0) * deg,
          f.number("C_nalpha", 0),
          f.number("T_f0", 3),
          f.number("T_V0", 6),
          f.number("T_p", 1.7),
          f.number("T_VL", 11),
          f.number("b1", .14),
          f.number("b2", .53),
          f.number("b5", 5),
          f.number("A1", .3),
          f.number("A2", .7),
          f.number("A5", 1),
          f.number("Cn1", 0),
          f.number("Cn2", 0),
          f.number("St_sh", .19),
          f.number("Cd0", 0),
          f.number("Cm0", 0),
          f.number("x_cp_bar", .2),
          f.number("UACutout", 45) * deg,
          f.number("filtCutOff", .5)};
    bool constant = true;
    for (const auto &x : af.coefficients)
        if (x.cl != af.coefficients.front().cl || x.cd != af.coefficients.front().cd ||
            x.cm != af.coefficients.front().cm)
            constant = false;
    enabled_ = !constant && f.flag("InclUAdata") && std::abs(p_.cnalpha) > 1e-12;
    if (c <= 0 || dt <= 0 || sound <= 0)
        throw std::runtime_error("Invalid UA dimensions");
}
std::pair<double, double> UnsteadyAirfoil::separation(double a, double slope) const {
    a = p_.alpha0 + std::remainder(a - p_.alpha0, 2 * pi);
    const auto c = af_->at(a);
    const double cn = c.cl * std::cos(a) + (c.cd - p_.cd0) * std::sin(a);
    const double ratio = small_single(slope) || small_single(a - p_.alpha0) || small_single(cn)
                             ? 0
                             : std::max(0., cn / (slope * (a - p_.alpha0)));
    const double f = std::min(1., sq(2 * std::sqrt(ratio) - 1));
    return {f, std::abs(f - 1) < 1e-14 ? 0 : (cn - slope * (a - p_.alpha0) * f) / (1 - f)};
}
double UnsteadyAirfoil::chord_separation(double a, double slope) const {
    a = std::remainder(a, 2 * pi);
    const double a0 = std::remainder(p_.alpha0, 2 * pi);
    if (small_single(a) || small_single(a - a0) || small_single(slope))
        return 1.44;
    const auto c = af_->at(a);
    const double cc = c.cl * std::sin(a) - (c.cd - p_.cd0) * std::cos(a);
    return std::min(1.44, sq(cc / (slope * std::remainder(a - a0, 2 * pi) * std::tan(a)) + .2));
}
UnsteadyAirfoil::Chain UnsteadyAirfoil::chain(double a, double speed) const {
    speed = std::max(.01, speed);
    a = std::remainder(a, 2 * pi);
    const auto &o = previous_;
    Chain k;
    const double M = speed / sound_;
    if (M >= 1)
        throw std::runtime_error("UA requires subsonic relative velocity");
    const double beta2 = 1 - M * M, beta = std::sqrt(beta2), oldalpha = first_ ? a : o.alpha;
    k.cna = p_.cnalpha / beta;
    k.ds = 2 * speed * dt_ / chord_;
    const double lp = std::exp(-2 * dt_ * std::max(1., speed) * p_.filter / chord_);
    k.alpha = lp * oldalpha + (1 - lp) * a;
    k.da0 = k.alpha - p_.alpha0;
    k.q = (k.alpha - oldalpha) / dt_ * chord_ / speed;
    const double oldq = first_ ? k.q : o.q, oldqf = first_ ? k.q : o.qf;
    k.qf = lp * oldqf + (1 - lp) * k.q;
    k.ka = k.qf * speed / chord_;
    k.kq = (k.q - oldq) / dt_;
    k.kqf = lp * (first_ ? 0 : o.kqf) + (1 - lp) * k.kq;
    const double ka = 1 / ((1 - M) + p_.cnalpha / 2 * M * M * beta * (p_.a1 * p_.b1 + p_.a2 * p_.b2));
    const double kq = 1 / ((1 - M) + p_.cnalpha * M * M * beta * (p_.a1 * p_.b1 + p_.a2 * p_.b2));
    const double TI = chord_ / sound_, Ta = TI * ka * .75, Tq = TI * kq * .75, Tf = p_.tf / sigma1_;
    k.kpa = decay(dt_, Ta, o.kpa, k.ka, first_ ? 0 : o.ka);
    k.cnanc = 4 * Ta * (k.ka - k.kpa) / M;
    k.kpq = decay(dt_, Tq, o.kpq, k.kqf, first_ ? 0 : o.kqf);
    k.cnqnc = -Tq * (k.kqf - k.kpq) / M;
    k.x1 = decay(k.ds * beta2 * p_.b1, 1, o.x1, p_.a1 * (k.alpha - oldalpha), 0);
    k.x2 = decay(k.ds * beta2 * p_.b2, 1, o.x2, p_.a2 * (k.alpha - oldalpha), 0);
    k.ae = k.da0 - k.x1 - k.x2;
    k.cncirc = k.cna * k.ae;
    k.k3q = decay(p_.b5 * beta2 * k.ds, 1, o.k3q, p_.a5 * (k.qf - oldqf), 0);
    k.cmqc = -p_.cnalpha * (k.qf - k.k3q) * chord_ / (16 * beta * speed);
    k.cnpot = k.cncirc + k.cnanc + k.cnqnc;
    const double kmq = 7 / (15 * (1 - M) + 1.5 * p_.cnalpha * p_.a5 * p_.b5 * beta * M * M);
    k.kppq = decay(dt_, kmq * kmq * TI, o.kppq, k.kqf, first_ ? 0 : o.kqf);
    k.cmqnc = -k.cnqnc / 4 - ka * ka * TI * (k.kqf - k.kppq) / (3 * M);
    k.ccpot = k.cncirc * std::tan(k.ae + p_.alpha0);
    k.dp = decay(k.ds, p_.tp, o.dp, k.cnpot, first_ ? k.cnpot : o.cnpot);
    k.cnprime = k.cnpot - k.dp;
    k.af = k.cnprime / k.cna + p_.alpha0;
    k.f = separation(k.af, k.cna).first;
    k.df = first_ ? 0 : decay(k.ds, Tf, o.df, k.f, o.f);
    k.fpp = k.f - k.df;
    k.fc = chord_separation(k.af, k.cna);
    k.dfc = first_ ? 0 : decay(k.ds, p_.tf, o.dfc, k.fc, o.fc);
    k.fcpp = k.fc - k.dfc;
    // Same domain as upstream: separation states must stay nonnegative.
    if (k.fpp < 0 || k.fcpp < 0)
        throw std::runtime_error("UA separation state became negative");
    k.cnfs =
        k.cnanc + k.cnqnc + k.cna * k.ae * k.fpp + separation(k.ae + p_.alpha0, k.cna).second * (1 - k.fpp);
    k.daf = first_ ? 0 : decay(k.ds, .1 * Tf, o.daf, k.af, o.af);
    k.cv = k.cncirc * (1 - sq(.5 + .5 * std::sqrt(k.fpp)));
    const double Tv = p_.tv / sigma3_;
    if (!first_)
        k.cnv = std::max(0., tau_ > p_.tvl && k.ka * k.da0 > 0 ? o.cnv * std::exp(-k.ds / Tv)
                                                               : decay(k.ds, Tv, o.cnv, k.cv, o.cv));
    return k;
}
Coefficients UnsteadyAirfoil::evaluate(double a, double speed) const {
    a = std::remainder(a, 2 * pi);
    const auto steady = af_->at(a);
    if (!enabled_ || first_)
        return steady;
    auto k = chain(a, speed);
    const double cn = k.cnfs + (tau_ > 0 ? k.cnv : 0);
    const double cc = k.ccpot * (std::sqrt(k.fcpp) - .2) + k.cnv * k.ae * (1 - tau_ / p_.tvl);
    Coefficients y{cn * std::cos(a) + cc * std::sin(a), cn * std::sin(a) - cc * std::cos(a) + p_.cd0,
                   af_->at(k.af - k.daf).cm + k.cmqc - k.cnanc / 4 + k.cmqnc -
                       p_.xcp * (1 - std::cos(pi * tau_ / p_.tvl)) * k.cnv};
    const double weight =
        (1 - blend(std::abs(a), p_.cutout - 5 * deg, p_.cutout)) * blend(std::abs(speed), .01, 1);
    return {weight * y.cl + (1 - weight) * steady.cl, weight * y.cd + (1 - weight) * steady.cd,
            weight * y.cm + (1 - weight) * steady.cm};
}
void UnsteadyAirfoil::advance(double a, double speed, std::size_t step) {
    if (!enabled_ || step == 0)
        return;
    const auto k = chain(a, speed);
    const bool les = k.cnprime > p_.cn1 || k.cnprime < p_.cn2, tes = k.fpp < previous_.fpp,
               vortex = tau_ <= p_.tvl && tau_ > 0;
    const double factor = k.ka * k.da0;
    sigma1_ = 1;
    if (tes) {
        if (factor < 0)
            sigma1_ = 2;
        else if (!les)
            sigma1_ = 1;
        else if (previous_.fpp <= .7)
            sigma1_ = 2;
        else
            sigma1_ = 1.75;
    } else {
        if (!les)
            sigma1_ = .5;
        if (vortex)
            sigma1_ = .25;
        if (factor > 0)
            sigma1_ = .75;
    }
    sigma3_ = 1;
    if (tau_ <= 2 * p_.tvl && tau_ >= p_.tvl) {
        sigma3_ = tes ? 3 : 4;
    } else if (vortex)
        sigma3_ = factor < 0 ? 2 : 1;
    else if (factor < 0)
        sigma3_ = 4;
    if (!tes && k.kqf * k.da0 < 0)
        sigma3_ = 1;
    if (tau_ > 0 || les)
        tau_ += k.ds;
    if (tau_ >= p_.tvl + 2 * (1 - k.fpp) / p_.st && tes)
        tau_ = 0;
    previous_ = k;
    first_ = false;
}
} // namespace turbine
