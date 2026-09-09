// BPM spectral shape functions from OpenFAST, Apache-2.0; see NOTICE.
#include "source_models.hpp"
#include <algorithm>
#include <cmath>
namespace aeroacoustics {
double g5comp(double thickness_ratio, double eta) {
    const double mu = thickness_ratio < .25    ? .1211
                      : thickness_ratio <= .62 ? -.2175 * thickness_ratio + .1755
                      : thickness_ratio < 1.15 ? -.0308 * thickness_ratio + .0596
                                               : .0242;
    double slope;
    if (thickness_ratio <= .02)
        slope = 0.;
    else if (thickness_ratio < .5)
        slope = 68.724 * thickness_ratio - 1.35;
    else if (thickness_ratio <= .62)
        slope = 308.475 * thickness_ratio - 121.23;
    else if (thickness_ratio <= 1.15)
        slope = 224.811 * thickness_ratio - 69.354;
    else if (thickness_ratio < 1.2)
        slope = 1583.28 * thickness_ratio - 1631.592;
    else
        slope = 268.344;
    slope = std::max(slope, 0.);
    const double eta0 = -std::sqrt(slope * slope * std::pow(mu, 4) / (6.25 + ((slope * slope) * mu) * mu));
    if (eta <= eta0) {
        const double offset = (2.5 * std::sqrt(1. - std::pow(eta0 / mu, 2)) - 2.5) - slope * eta0;
        return slope * eta + offset;
    }
    if (eta <= 0.)
        return 2.5 * std::sqrt(1. - std::pow(eta / mu, 2)) - 2.5;
    if (eta <= .03615995)
        return std::sqrt(1.5625 - 1194.99 * std::pow(eta, 2)) - 1.25;
    return -155.543 * eta + 4.375;
}
double amin(double argument) {
    const double x = std::abs(argument);
    if (x <= .204)
        return std::sqrt(67.552 - 886.788 * std::pow(x, 2)) - 8.219;
    if (x <= .244)
        return -32.665 * x + 3.981;
    return ((-142.795 * std::pow(x, 3) + 103.656 * std::pow(x, 2)) - 57.757 * x) + 6.006;
}
double amax(double argument) {
    const double x = std::abs(argument);
    if (x <= .13)
        return std::sqrt(67.552 - 886.788 * std::pow(x, 2)) - 8.219;
    if (x <= .321)
        return -15.901 * x + 1.098;
    return ((-4.669 * std::pow(x, 3) + 3.491 * std::pow(x, 2)) - 16.699 * x) + 1.149;
}
double bmin(double argument) {
    const double x = std::abs(argument);
    if (x <= .13)
        return std::sqrt(16.888 - 886.788 * std::pow(x, 2)) - 4.109;
    if (x <= .145)
        return -83.607 * x + 8.138;
    return ((-817.81 * std::pow(x, 3) + 355.21 * std::pow(x, 2)) - 135.024 * x) + 10.619;
}
double bmax(double argument) {
    const double x = std::abs(argument);
    if (x <= .1)
        return std::sqrt(16.888 - 886.788 * std::pow(x, 2)) - 4.109;
    if (x <= .187)
        return -31.313 * x + 1.854;
    return ((-80.541 * std::pow(x, 3) + 44.174 * std::pow(x, 2)) - 39.381 * x) + 2.344;
}
double a0comp(double reynolds) {
    if (reynolds < 95200.)
        return .57;
    if (reynolds < 857000.)
        return -9.57e-13 * std::pow(reynolds - 857000., 2) + 1.13;
    return 1.13;
}
} // namespace aeroacoustics