#pragma once
#include "aeroacoustics.hpp"
namespace aeroacoustics {
double log10aa(double x);
Spectrum lblvs(double alpstar, double c, double u, double theta, double phi, double l, double r, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval);
std::tuple<Spectrum, Spectrum, Spectrum> tblte(double alpstar, double c, double u, double theta, double phi, double l, double r, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval);
Spectrum tipnois(double alphtip, double alprat2, double c, double u, double theta, double phi, double r, const Parameters& p);
Spectrum inflownoise(double alphanoise, double chord, double u, double theta, double phi, double d, double robs, double tinoise, const Parameters& p);
Spectrum blunt(double alpstar, double c, double u, double theta, double phi, double l, double r, double h, double psi, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval);
double g5comp(double hdstar, double eta);
double amin(double a);
double amax(double a);
double bmin(double b);
double bmax(double b);
double a0comp(double rc);
std::tuple<double, double, double> thick(double c, double rc, double alpstar, const Parameters& p, double stallval);
double directh_te(double m, double theta, double phi);
double directh_le(double m, double theta, double phi);
double directl(double m, double theta, double phi);
Spectrum simple_guidati(double u, double chord, double thick_10p, double thick_1p, const Parameters& p);
}
