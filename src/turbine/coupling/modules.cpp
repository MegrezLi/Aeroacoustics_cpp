#include "turbine/modules.hpp"
#include "turbine/input.hpp"
#include <algorithm>
namespace turbine {
const std::array<ModuleCapability, 6> &module_capabilities() {
    static const std::array<ModuleCapability, 6> capabilities{
        {{ModuleId::steady_wind, "steady_wind", "global position", "wind velocity [m/s]", "none"},
         {ModuleId::fixed_base_blades, "fixed_base_blades", "modal q/qd and frozen aerodynamic loads",
          "modal acceleration; mesh motion", "q/qd and integrator acceleration history"},
         {ModuleId::steady_bem, "steady_bem", "predicted motion, wind, polars",
          "station loads, speed, incidence", "previous phi root"},
         {ModuleId::minnema_pierce_ua, "minnema_pierce_ua", "previous aerodynamic station state",
          "unsteady airfoil coefficients", "UA_Mod=3 discrete states"},
         {ModuleId::generalized_alpha, "generalized_alpha", "DOF layout and acceleration operator",
          "corrected q/qd", "physical and algorithmic acceleration"},
         {ModuleId::acoustics, "aeroacoustics", "predicted station state, BL, observers",
          "band sound pressure level", "TI buffers and sampling time"}}};
    return capabilities;
}
DofLayout fixed_base_dof_layout(std::array<double, 3> scales) {
    std::vector<DofDescriptor> dofs;
    std::vector<DofBlock> blocks;
    const std::array<const char *, 3> names{"flap1", "flap2", "edge"};
    for (std::size_t b = 0; b < FixedBaseBladeBackend::blades; ++b) {
        blocks.push_back({"blade" + std::to_string(b + 1), dofs.size(), names.size()});
        for (std::size_t j = 0; j < names.size(); ++j)
            dofs.push_back({"blade" + std::to_string(b + 1) + "_" + names[j], "m", scales[j],
                            std::to_string(b + 1) + "_" + names[j]});
    }
    return {std::move(dofs), std::move(blocks)};
}
ModuleConfiguration configure_modules(const Case &c) {
    const auto require = [](bool ok, const std::string &message) {
        if (!ok)
            throw std::runtime_error(message);
    };
    const auto &primary = c.primary, &structure = c.structure, &aero = c.aero;
    require(primary.integer("CompElast") == 1 && primary.integer("CompInflow") == 1 &&
                primary.integer("CompAero") == 2,
            "Unsupported module configuration: expected fixed-base blades + steady wind + AeroDyn BEM");
    for (auto key :
         {"CompServo", "CompSeaSt", "CompHydro", "CompSub", "CompMooring", "CompIce", "CompSoil", "MHK"})
        require(primary.integer(key) == 0, std::string("Unsupported enabled module: ") + key);
    require(primary.integer("NRotors") == 1 && !primary.flag("MirrorRotor"),
            "Only one normal rotor is supported");
    require(structure.integer("NumBl") == FixedBaseBladeBackend::blades,
            "Only three blades are supported by fixed_base_blades");
    for (auto key : {"FlapDOF1", "FlapDOF2", "EdgeDOF"})
        require(structure.flag(key), std::string("Expected active blade DOF: ") + key);
    for (auto key :
         {"PitchDOF", "TeetDOF", "DrTrDOF", "GenDOF", "YawDOF", "TwFADOF1", "TwFADOF2", "TwSSDOF1",
          "TwSSDOF2", "PtfmSgDOF", "PtfmSwDOF", "PtfmHvDOF", "PtfmRDOF", "PtfmPDOF", "PtfmYDOF", "Furling"})
        require(!structure.flag(key), std::string("Unsupported active DOF: ") + key +
                                          " (requires another structural backend and coupled DOF layout)");
    require(c.inflow.integer("WindType") == 1, "Only steady InflowWind WindType=1 is supported");
    require(aero.integer("Wake_Mod") == 1 && aero.integer("BEM_Mod") == 1 && aero.integer("DBEMT_Mod") == 0,
            "Unsupported wake model");
    require(aero.integer("UA_Mod") == 3 && aero.flag("FLookup") && aero.integer("AFTabMod") == 1,
            "Unsupported unsteady airfoil model");
    require(primary.integer("ModCoupling") == 3 && primary.integer("NumCrctn") == 0,
            "Expected ModCoupling=3 and NumCrctn=0");
    return {{{ModuleId::steady_wind, ModuleId::fixed_base_blades, ModuleId::steady_bem,
              ModuleId::minnema_pierce_ua, ModuleId::generalized_alpha, ModuleId::acoustics}},
            FixedBaseBladeBackend::blades};
}
DofLayout structural_dof_layout(const ModuleConfiguration &config, std::array<double, 3> scales) {
    for (const auto &capability : module_capabilities())
        if (std::count(config.modules.begin(), config.modules.end(), capability.id) != 1)
            throw std::invalid_argument("Unsupported or duplicate module combination");
    if (config.blades != FixedBaseBladeBackend::blades)
        throw std::invalid_argument("Selected structural backend requires three blades");
    return fixed_base_dof_layout(scales);
}
} // namespace turbine
