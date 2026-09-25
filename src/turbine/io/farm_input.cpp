#include "turbine/engineering.hpp"
#include "turbine/farm.hpp"
namespace turbine {
FarmOptions FarmOptions::read(const std::filesystem::path &p) {
    InputFile f(p);
    FarmOptions o;
    o.provenance = f.value("Provenance");
    o.duration = f.number("Duration");
    o.statistics_start = f.number("StatisticsStart");
    o.wakes = f.flag("Wakes");
    o.expansion = f.number("WakeExpansion");
    o.max_ct = f.number("MaxCt");
    int memory = f.integer("MaxValues");
    if (memory <= 0)
        throw std::invalid_argument("Invalid farm MaxValues");
    o.max_values = std::size_t(memory);
    o.observers = aeroacoustics::read_observers(f.file("ObserverFile").string());
    if (f.value("Propagation") != "none")
        o.propagation = read_propagation(f.file("Propagation"));
    int nt = f.integer("NumTurbines"), ns = f.integer("NumSources");
    if (nt < 1 || nt > 1000 || ns < 0 || ns > 10000)
        throw std::invalid_argument("Invalid farm counts");
    for (const auto &path : f.files_after("TurbineFiles", nt)) {
        InputFile u(path);
        FarmUnit unit;
        unit.name = u.value("Name");
        unit.case_file = u.file("Case");
        unit.origin = {u.number("X"), u.number("Y"), 0};
        if (u.value("Controller") != "none")
            unit.controller = ControlConfig::read(u.file("Controller"));
        if (u.value("Tower") != "none")
            unit.tower = TowerInfluence::read(u.file("Tower"));
        if (u.value("Surfaces") != "none")
            unit.surfaces = SurfaceSet::read(u.file("Surfaces"));
        o.turbines.push_back(std::move(unit));
    }
    if (ns)
        for (const auto &path : f.files_after("SourceFiles", ns))
            o.sources.push_back(StationarySource::read(path));
    return o;
}
} // namespace turbine
