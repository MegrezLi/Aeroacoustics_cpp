#pragma once
#include <array>
#include <cstddef>
namespace aeroacoustics {
enum class Mechanism : std::size_t {
    laminar,
    trailing_pressure,
    trailing_suction,
    trailing_separation,
    bluntness,
    tip,
    inflow,
    count
};
inline constexpr std::size_t mechanism_count = static_cast<std::size_t>(Mechanism::count);
constexpr std::size_t index(Mechanism m) noexcept { return static_cast<std::size_t>(m); }
struct MechanismInfo {
    Mechanism id;
    const char *name;
    const char *unit;
    bool contributes_to_total;
};
// Stable legacy channel order. Models replace contributions within these mechanisms.
inline constexpr std::array<MechanismInfo, mechanism_count> mechanism_registry{
    {{Mechanism::laminar, "LBL", "dB", true},
     {Mechanism::trailing_pressure, "TBL_pressure", "dB", true},
     {Mechanism::trailing_suction, "TBL_suction", "dB", true},
     {Mechanism::trailing_separation, "TBL_separation", "dB", true},
     {Mechanism::bluntness, "bluntness", "dB", true},
     {Mechanism::tip, "tip", "dB", true},
     {Mechanism::inflow, "inflow", "dB", true}}};
constexpr bool valid_registry() {
    for (std::size_t i = 0; i < mechanism_registry.size(); ++i)
        if (index(mechanism_registry[i].id) != i)
            return false;
    return true;
}
static_assert(valid_registry(), "Mechanism descriptors must follow indexed storage order");
} // namespace aeroacoustics
