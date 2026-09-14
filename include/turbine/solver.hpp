#pragma once
#include "turbine/rotor.hpp"
namespace turbine {
class Solver {
  public:
    explicit Solver(const Case &);
    explicit Solver(TurbineModel);
    // Copies share frozen model data; all evolving states and workspaces are independent.
    const Rotor &rotor() const noexcept { return rotor_; }
    const RotorState &state() const noexcept { return state_; }
    const RotorOutput &aerodynamic() const noexcept { return aerodynamic_; }
    const std::array<Vec3, 3> &acceleration() const noexcept { return acceleration_; }
    double time() const noexcept { return time_; }
    std::size_t step_number() const noexcept { return step_number_; }
    bool failed() const noexcept { return failed_; }
    void step();
    void reset();
    class Checkpoint {
        friend class Solver;
        std::shared_ptr<const Solver> saved_;
        explicit Checkpoint(const Solver &s) : saved_(std::make_shared<const Solver>(s)) {}

      public:
        Checkpoint(const Checkpoint &) = default;
        Checkpoint &operator=(const Checkpoint &) = default;
    };
    Checkpoint checkpoint() const;
    void restore(const Checkpoint &);

  private:
    Rotor rotor_;
    RotorState state_{};
    RotorOutput aerodynamic_;
    std::array<Vec3, 3> acceleration_{}, algorithmic_{};
    double time_ = 0;
    std::size_t step_number_ = 0;
    bool failed_ = false;
    RotorWorkspace rotor_workspace_;
    LoadWorkspace load_workspace_;
    double dt_, alpha_m_, alpha_f_, beta_, gamma_, beta_prime_, gamma_prime_;
    void advance();
};
} // namespace turbine
