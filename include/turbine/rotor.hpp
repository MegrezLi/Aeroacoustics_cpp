#pragma once
#include "turbine/bem.hpp"
#include "turbine/mesh.hpp"
#include "turbine/model.hpp"
#include "turbine/unsteady.hpp"
#include "turbine/wind.hpp"
namespace turbine {
using RotorState = std::array<ModalState, FixedBaseBladeBackend::blades>;
struct AeroStation {
    Motion motion;
    Vec3 wind;
    Matrix3 annulus;
    BEMInput bem;
    double root_phi = 0, phi = 0, alpha = 0, speed = 0, axial = 0, tangential = 0;
    Coefficients coefficients;
    PointLoad load;
};
struct RotorOutput {
    std::array<std::vector<AeroStation>, FixedBaseBladeBackend::blades> blades;
    Vec3 average_velocity{};
    double skew = 0;
};
// Owned by the caller/Solver, never shared as hidden mutable Rotor state.
struct RotorWorkspace {
    std::vector<Motion> structural, aerodynamic;
};
struct LoadWorkspace {
    std::vector<Vec3> source, destination;
    std::vector<PointLoad> distributed, points, result;
};
class Rotor {
    TurbineModel model_;
    BladeStructure structure_;

  public:
    explicit Rotor(const Case &);
    explicit Rotor(TurbineModel);
    const TurbineModel &model() const noexcept { return model_; }
    const BladeStructure &structure() const noexcept { return structure_; }
    void set_operation(const RotorKinematics &op) { structure_.set_operation(op); }
    void set_wind(std::shared_ptr<const WindField> wind) { wind_ = std::move(wind); }
    Vec3 wind_at(double time, const Vec3 &p) const {
        return wind_ ? wind_->at(time, p) : model_.data().wind.at(p);
    }
    double aerodynamic_torque(const RotorOutput &) const;
    RotorOutput evaluate(double time, const RotorState &) const;
    void evaluate_into(double time, const RotorState &, RotorOutput &, RotorWorkspace &) const;
    // Reuses motions for both load mapping and acceleration in one state evaluation.
    Vec3 structural_acceleration(std::size_t blade, const Matrix3 &basis, const ModalState &,
                                 const RotorOutput &, std::vector<Motion> &, LoadWorkspace &) const;
    void advance_airfoils(const RotorOutput &, std::size_t step);
    std::array<std::vector<PointLoad>, 3> structural_loads(double time, const RotorState &,
                                                           const RotorOutput &) const;
    std::vector<PointLoad> structural_loads_for_blade(std::size_t blade, double time, const ModalState &,
                                                      const RotorOutput &) const;

  private:
    void structural_loads_into(std::size_t blade, const RotorOutput &, const std::vector<Motion> &,
                               LoadWorkspace &) const;
    BEMOptions options_;
    std::shared_ptr<const WindField> wind_;
    struct SkewOptions {
        bool redistribute;
        double factor;
    } skew_;
    std::vector<double> tip_constant_, hub_constant_;
    std::vector<MotionMap> motion_maps_;
    std::vector<LoadMap> load_maps_;
    std::array<std::vector<UnsteadyAirfoil>, FixedBaseBladeBackend::blades> airfoils_;
    std::array<std::vector<double>, FixedBaseBladeBackend::blades> previous_phi_;
};
} // namespace turbine
