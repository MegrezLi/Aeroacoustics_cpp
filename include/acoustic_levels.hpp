#pragma once
#include "aeroacoustics.hpp"
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace aeroacoustics {
inline void validate_level(double level) {
    if (std::isnan(level) || level == std::numeric_limits<double>::infinity())
        throw std::runtime_error("Invalid acoustic level (NaN or +Inf)");
}
inline double relative_power(double level) {
    validate_level(level);
    if (level == -std::numeric_limits<double>::infinity())
        return 0.;
    const double power = std::pow(10., level / 10.);
    if (!std::isfinite(power))
        throw std::runtime_error("Acoustic power overflow");
    return power;
}
inline double output_decibels(double power) {
    if (!std::isfinite(power) || power < 0)
        throw std::runtime_error("Invalid accumulated acoustic power");
    // Legacy Fortran-compatible zero placeholder; accompanied by a zero mask.
    return power == 0 ? 0 : 10 * std::log10(power);
}
inline std::string spectrum_context(const Parameters &p, std::size_t mechanism, std::size_t frequency) {
    static const char *names[] = {"LBL",       "TBL_pressure", "TBL_suction", "TBL_separation",
                                  "bluntness", "tip",          "inflow"};
    std::ostringstream message;
    message << " mechanism=" << names[mechanism] << " frequency_Hz=" << p.freqlist[frequency]
            << " TBLTEMod=" << p.tbltemod << " TIMod=" << p.timod;
    return message.str();
}
} // namespace aeroacoustics
