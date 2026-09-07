#pragma once
#include "turbine/input.hpp"
namespace turbine {
// Direct port of the UA_Mod=3 Kelvin chain and discrete states. evaluate()
// leaves state untouched; advance() is called once per aerodynamic time step.
class UnsteadyAirfoil {
  public:
    UnsteadyAirfoil(const Airfoil &, double chord, double dt, double sound_speed);
    Coefficients evaluate(double alpha, double speed) const;
    void advance(double alpha, double speed, std::size_t step);
    bool active() const { return enabled_; }

  private:
    const Airfoil *af_;
    double chord_, dt_, sound_;
    bool enabled_, first_ = true;
    struct Parameters {
        double alpha0, cnalpha, tf, tv, tp, tvl, b1, b2, b5, a1, a2, a5, cn1, cn2, st, cd0, cm0, xcp, cutout,
            filter;
    } p_;
    struct Chain {
        double alpha = 0, q = 0, qf = 0, ka = 0, kq = 0, kqf = 0, x1 = 0, x2 = 0, kpa = 0, kpq = 0, k3q = 0,
               kppq = 0;
        double dp = 0, cnpot = 0, cnprime = 0, f = 0, fpp = 0, df = 0, fc = 0, fcpp = 0, dfc = 0, daf = 0,
               af = 0, cnv = 0, cv = 0;
        double ds = 0, cncirc = 0, cnanc = 0, cnqnc = 0, cna = 0, ae = 0, cnfs = 0, ccpot = 0, cmqc = 0,
               cmqnc = 0, da0 = 0;
    } previous_;
    double sigma1_ = 1, sigma3_ = 1, tau_ = 0;
    Chain chain(double alpha, double speed) const;
    std::pair<double, double> separation(double alpha, double slope) const;
    double chord_separation(double alpha, double slope) const;
};
} // namespace turbine
