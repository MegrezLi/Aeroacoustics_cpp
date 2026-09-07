#pragma once
#include "turbine/input.hpp"
namespace turbine {
struct BEMOptions {
    int blades = 3;
    bool tip_loss = true, hub_loss = true, tangential = true, axial_drag = false, tangential_drag = false;
    double tolerance = 5e-10;
    int max_iterations = 100;
};
struct BEMInput {
    double radius, chord, twist, vx, vy, hub_loss_constant, tip_loss_constant;
};
struct Induction {
    double axial = 0, tangential = 0, residual = 0, loss = 1, phi = 0, k = 0, kp = 0;
    bool valid = true;
};
double prandtl_loss(const BEMOptions &, const BEMInput &, double phi);
Induction induction_factors(const BEMOptions &, const BEMInput &, double phi, double cn, double ct,
                            double loss);
Induction bem_residual(const BEMOptions &, const BEMInput &, const Airfoil &, double phi);
Induction solve_bem(const BEMOptions &, const BEMInput &, const Airfoil &, double previous_phi);
double skew_axial(double axial, double chi0, double radius_ratio, double azimuth,
                  double factor = 15 * pi / 32);
} // namespace turbine
