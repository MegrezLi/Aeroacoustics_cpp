// Directivity correlations from OpenFAST, Apache-2.0; see NOTICE.
#include "kernels.hpp"
#include <cmath>

namespace aeroacoustics {
namespace {
// Keep the coefficient of the reference empirical model, rather than pi/180.
constexpr double directivity_degrees = .017453;
} // namespace
double directh_te(double mach, double theta, double phi) {
    const double angle = theta * directivity_degrees, azimuth = phi * directivity_degrees;
    const double convection_mach = .8 * mach;
    return 2. * std::pow(std::sin(angle / 2.), 2) * std::pow(std::sin(azimuth), 2) /
           ((1. + mach * std::cos(angle)) * std::pow(1. + (mach - convection_mach) * std::cos(angle), 2));
}
double directh_le(double mach, double theta, double phi) {
    const double angle = theta * directivity_degrees, azimuth = phi * directivity_degrees;
    return 2. * std::pow(std::cos(angle / 2.), 2) * std::pow(std::sin(azimuth), 2) /
           std::pow(1. + mach * std::cos(angle), 3);
}
double directl(double mach, double theta, double phi) {
    const double angle = theta * directivity_degrees, azimuth = phi * directivity_degrees;
    return std::pow(std::sin(angle) * std::sin(azimuth), 2) / std::pow(1. + mach * std::cos(angle), 4);
}
} // namespace aeroacoustics
