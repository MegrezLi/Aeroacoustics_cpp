#include "turbine/engineering.hpp"
namespace turbine {
aeroacoustics::PropagationOptions read_propagation(const std::filesystem::path &path) {
    InputFile f(path);
    aeroacoustics::PropagationOptions p;
    p.absorption = f.flag("AirAbsorption");
    p.atmosphere = {f.number("TemperatureC") + 273.15, f.number("RelativeHumidity"), f.number("PressurePa")};
    p.sound_speed = f.number("SoundSpeed");
    p.ground_z = f.number("GroundZ");
    const auto ground = f.value("GroundModel");
    if (ground == "none")
        p.ground = aeroacoustics::GroundModel::none;
    else if (ground == "rigid")
        p.ground = aeroacoustics::GroundModel::rigid;
    else if (ground == "impedance")
        p.ground = aeroacoustics::GroundModel::impedance;
    else
        throw std::invalid_argument("GroundModel must be none, rigid or impedance");
    p.normalized_impedance = {f.number("ImpedanceReal"), f.number("ImpedanceImag")};
    const int count = f.integer("NumScreens");
    if (count < 0 || count > 10000)
        throw std::invalid_argument("Invalid screen count");
    if (count)
        for (const auto &v : f.table_after("NumScreens", count, 5))
            p.screens.push_back({v[0], v[1], v[2], v[3], v[4]});
    p.validate();
    return p;
}
} // namespace turbine
