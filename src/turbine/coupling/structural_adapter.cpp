#include "turbine/structural_adapter.hpp"
namespace turbine {
FixedBaseAcceleration::FixedBaseAcceleration(const Rotor &rotor, const RotorOutput &aero, double time,
                                             RotorWorkspace &motions, LoadWorkspace &loads)
    : rotor_(rotor), aerodynamic_(aero), motions_(motions), loads_(loads) {
    for (std::size_t b = 0; b < basis_.size(); ++b)
        basis_[b] = rotor_.structure().blade_basis(time, b);
}
void FixedBaseAcceleration::evaluate(std::size_t b, StateView state, double *out) const {
    if (b >= basis_.size() || state.size != FixedBaseBladeBackend::modes_per_blade)
        throw std::invalid_argument("Fixed-base acceleration layout mismatch");
    ModalState modal;
    std::copy_n(state.q, modal.q.size(), modal.q.begin());
    std::copy_n(state.qd, modal.qd.size(), modal.qd.begin());
    const auto f =
        rotor_.structural_acceleration(b, basis_[b], modal, aerodynamic_, motions_.structural, loads_);
    std::copy(f.begin(), f.end(), out);
}
RotorState rotor_state(const SecondOrderState &state) {
    constexpr auto modes = FixedBaseBladeBackend::modes_per_blade;
    if (state.q.size() != FixedBaseBladeBackend::blades * modes || state.qd.size() != state.q.size())
        throw std::invalid_argument("Fixed-base state layout mismatch");
    RotorState result;
    for (std::size_t b = 0; b < result.size(); ++b) {
        std::copy_n(state.q.data() + b * modes, modes, result[b].q.begin());
        std::copy_n(state.qd.data() + b * modes, modes, result[b].qd.begin());
    }
    return result;
}
std::array<Vec3, FixedBaseBladeBackend::blades> rotor_acceleration(const std::vector<double> &a) {
    constexpr auto modes = FixedBaseBladeBackend::modes_per_blade;
    if (a.size() != FixedBaseBladeBackend::blades * modes)
        throw std::invalid_argument("Fixed-base acceleration shape mismatch");
    std::array<Vec3, FixedBaseBladeBackend::blades> result;
    for (std::size_t b = 0; b < result.size(); ++b)
        std::copy_n(a.data() + b * modes, modes, result[b].begin());
    return result;
}
} // namespace turbine
