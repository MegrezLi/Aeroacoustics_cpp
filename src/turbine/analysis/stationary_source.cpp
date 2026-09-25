#include "acoustic_levels.hpp"
#include "turbine/farm.hpp"
#include <algorithm>
namespace turbine {
void StationarySource::validate() const {
    if (name.empty() || provenance.empty() || (kind != "bands" && kind != "tones") || frequencies.empty() ||
        frequencies.size() != levels_db.size() || !std::isfinite(directivity) || std::abs(directivity) > 1)
        throw std::invalid_argument("Invalid stationary source metadata/spectrum");
    for (double x : position)
        if (!std::isfinite(x))
            throw std::invalid_argument("Invalid source position");
    for (double x : axis)
        if (!std::isfinite(x))
            throw std::invalid_argument("Invalid source axis");
    if (std::abs(norm(axis) - 1) > 1e-8)
        throw std::invalid_argument("Source axis must have unit length");
    for (std::size_t i = 0; i < frequencies.size(); ++i)
        if (!std::isfinite(frequencies[i]) || frequencies[i] <= 0 || !std::isfinite(levels_db[i]) ||
            (i && frequencies[i] <= frequencies[i - 1]))
            throw std::invalid_argument("Invalid source frequency/level");
}
StationarySource StationarySource::read(const std::filesystem::path &p) {
    InputFile f(p);
    StationarySource s;
    s.name = f.value("Name");
    s.provenance = f.value("Provenance");
    s.kind = f.value("Kind");
    s.position = {f.number("X"), f.number("Y"), f.number("Z")};
    s.axis = {f.number("AxisX"), f.number("AxisY"), f.number("AxisZ")};
    s.directivity = f.number("Directivity");
    const int n = f.integer("NumFrequencies");
    if (n < 1 || n > 10000)
        throw std::invalid_argument("Source frequency count");
    for (const auto &r : f.table_after("Spectrum", n, 2)) {
        s.frequencies.push_back(r[0]);
        s.levels_db.push_back(r[1]);
    }
    s.validate();
    return s;
}
std::vector<double>
StationarySource::receiver_power(const Vec3 &receiver, const aeroacoustics::FrequencyBands &bands, double rho,
                                 double c,
                                 const std::optional<aeroacoustics::PropagationOptions> &prop) const {
    validate();
    for (double x : receiver)
        if (!std::isfinite(x))
            throw std::invalid_argument("Invalid source receiver");
    if (!std::isfinite(rho) || !std::isfinite(c) || rho <= 0 || c <= 0)
        throw std::invalid_argument("Invalid acoustic medium");
    if (prop) {
        prop->validate();
        if (prop->ground != aeroacoustics::GroundModel::none || !prop->screens.empty())
            throw std::invalid_argument("Stationary source supports free field/air absorption only");
    }
    const Vec3 delta = receiver - position;
    const double distance = norm(delta);
    if (!std::isfinite(distance) || distance <= 0)
        throw std::invalid_argument("Receiver coincides with stationary source");
    const double q = std::max(0., 1 + directivity * dot(delta / distance, axis));
    std::vector<double> out(bands.size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        const double f = frequencies[i];
        std::size_t bin = bands.size();
        for (std::size_t j = 0; j < bands.size(); ++j) {
            const auto &b = bands.values()[j];
            if ((kind == "bands" && std::abs(f - b.center_hz) < 1e-8) ||
                (kind == "tones" && f >= b.lower_hz &&
                 (f < b.upper_hz || (j + 1 == bands.size() && f <= b.upper_hz)))) {
                bin = j;
                break;
            }
        }
        if (bin == bands.size())
            throw std::invalid_argument("Stationary source frequency outside/unequal to requested bands");
        const double loss = prop && prop->absorption ? prop->atmosphere.absorption_db_per_m(f) * distance : 0;
        // p^2 = rho*c * W*Q/(4*pi*r^2); Lw is unweighted dB re 1 pW.
        const double lp = levels_db[i] + 10 * std::log10(rho * c * 1e-12 / (4 * pi * 4e-10)) -
                          20 * std::log10(distance) + aeroacoustics::a_weighting({f})[0] - loss;
        const double power = q == 0 ? 0 : q * aeroacoustics::relative_power(lp);
        out[bin] += power;
        if (!std::isfinite(out[bin]))
            throw std::overflow_error("Stationary source energy overflow");
    }
    return out;
}
} // namespace turbine
