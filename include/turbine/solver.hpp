#pragma once
#include "turbine/rotor.hpp"
namespace turbine {
class Solver {
  public:
    explicit Solver(const Case &);
    Rotor rotor;
    RotorState state{};
    RotorOutput aerodynamic;
    std::array<Vec3, 3> acceleration{}, algorithmic{};
    double time = 0;
    std::size_t step_number = 0;
    void step();

  private:
    double dt_, alpha_m_, alpha_f_, beta_, gamma_, beta_prime_, gamma_prime_;
};
} // namespace turbine
