#include "turbine/metrics.hpp"
namespace turbine {
void MetricsOptions::validate() const {
    for (double x : {end, receiver_dt, am_window, am_min_hz, am_max_hz, wind_bin_width})
        if (!std::isfinite(x) || x <= 0)
            throw std::invalid_argument("Invalid positive metrics parameter");
    if (!std::isfinite(start) || start < 0 || end <= start || am_max_hz < am_min_hz ||
        am_max_hz >= .5 / receiver_dt || am_window < 2 / am_min_hz || !max_values)
        throw std::invalid_argument("Invalid metrics window/Nyquist/memory limit");
    if (background_laeq_db && !std::isfinite(*background_laeq_db))
        throw std::invalid_argument("Invalid background level");
    if ((end - start) / receiver_dt > 10000000 || am_window / receiver_dt > 10000000)
        throw std::invalid_argument("Excessive receiver grid");
    if (std::abs(am_window / receiver_dt - std::round(am_window / receiver_dt)) > 1e-8)
        throw std::invalid_argument("AMWindow must be an integer number of ReceiverDT samples");
    for (const auto &p : observers)
        for (double v : p)
            if (!std::isfinite(v))
                throw std::invalid_argument("Invalid receiver map");
    for (const auto &t : tones) {
        if (t.name.empty() || t.provenance.empty() || !std::isfinite(t.frequency_hz) ||
            !std::isfinite(t.rotor_order) || t.frequency_hz < 0 || t.rotor_order < 0 ||
            (t.frequency_hz == 0 && t.rotor_order == 0) || !std::isfinite(t.level_db) ||
            !std::isfinite(t.reference_distance) || t.reference_distance <= 0)
            throw std::invalid_argument("Invalid tone input/provenance");
    }
}
MetricsOptions MetricsOptions::read(const std::filesystem::path &path) {
    InputFile f(path);
    MetricsOptions o;
    o.start = f.number("Start");
    o.end = f.number("End");
    o.receiver_dt = f.number("ReceiverDT");
    o.am_window = f.number("AMWindow");
    o.am_min_hz = f.number("AMMinHz");
    o.am_max_hz = f.number("AMMaxHz");
    o.wind_bin_width = f.number("WindBinWidth");
    o.retarded_time = f.flag("RetardedTime");
    o.apparent_power = f.flag("ApparentPower");
    const int memory = f.integer("MaxValues");
    if (memory <= 0)
        throw std::invalid_argument("Invalid MaxValues");
    o.max_values = std::size_t(memory);
    if (f.value("BackgroundLAeq") != "none")
        o.background_laeq_db = f.number("BackgroundLAeq");
    if (f.value("ObserverFile") != "none")
        o.observers = aeroacoustics::read_observers(f.file("ObserverFile").string());
    int nt = f.integer("NumTones");
    if (nt < 0 || nt > 10000)
        throw std::invalid_argument("Invalid tone count");
    if (nt)
        for (const auto &file : f.files_after("ToneFiles", std::size_t(nt))) {
            InputFile tone(file);
            ToneSource t;
            t.name = tone.value("Name");
            t.provenance = tone.value("Provenance");
            int b = tone.integer("Blade"), n = tone.integer("Node");
            if (b < 1 || n < 1)
                throw std::invalid_argument("Tone location must be one based");
            t.blade = std::size_t(b - 1);
            t.node = std::size_t(n - 1);
            t.frequency_hz = tone.number("FrequencyHz");
            t.rotor_order = tone.number("RotorOrder");
            t.level_db = tone.number("ReferenceSPL");
            t.reference_distance = tone.number("ReferenceDistance");
            o.tones.push_back(t);
        }
    o.validate();
    return o;
}
} // namespace turbine
