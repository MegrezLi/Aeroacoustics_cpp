#pragma once
#include "aeroacoustics.hpp"
namespace aeroacoustics {
// Howe-Chase, Mayer et al. (2019), equations 16-19. q = omega*delta/Uc.
double howe_chase_shape(double q, double half_height, double wavelength, double delta);
// One-sided PSD in Pa^2/Hz; theta/phi follow the existing trailing-edge Geometry convention.
double howe_chase_psd(double frequency_hz, double edge_speed, double delta, double span, const Geometry &,
                      const Parameters &);
struct TrailingEdgeOptions {
    std::string provenance;
    int model = 1;
    TnoEdgeVelocity edge_velocity = TnoEdgeVelocity::reference;
    HoweOptions howe;
    void validate() const;
    void apply(Parameters &) const;
    static TrailingEdgeOptions read(const std::string &path);
};
} // namespace aeroacoustics
