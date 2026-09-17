#pragma once
#include "turbine/integrator.hpp"
#include <array>
namespace turbine {
struct Case;
// Compile-time dimensions belong to this backend, not the numerical integrator.
struct FixedBaseBladeBackend {
    static constexpr std::size_t blades = 3, modes_per_blade = 3;
};
enum class ModuleId {
    steady_wind,
    fixed_base_blades,
    steady_bem,
    minnema_pierce_ua,
    generalized_alpha,
    acoustics
};
struct ModuleCapability {
    ModuleId id;
    const char *name, *inputs, *outputs, *history;
};
const std::array<ModuleCapability, 6> &module_capabilities();
struct ModuleConfiguration {
    std::array<ModuleId, 6> modules;
    std::size_t blades;
    const char *profile = "fixed_base_three_blade_three_mode";
};
ModuleConfiguration configure_modules(const Case &);
DofLayout fixed_base_dof_layout(std::array<double, 3> acceleration_scales = {1., 1., 1.});
DofLayout structural_dof_layout(const ModuleConfiguration &,
                                std::array<double, 3> acceleration_scales = {1., 1., 1.});
} // namespace turbine
