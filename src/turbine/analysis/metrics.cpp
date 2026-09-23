#include "turbine/metrics.hpp"
#include "acoustic_levels.hpp"
#include "checked_output.hpp"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <sstream>
namespace turbine {
namespace {
double distance(const Vec3 &a, const Vec3 &b) {
    return std::hypot(std::hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}
double level(double p) {
    if (!std::isfinite(p) || p < 0)
        throw std::invalid_argument("Invalid metric energy");
    return p == 0 ? -INFINITY : 10 * std::log10(p);
}
std::string quote(const std::string &s) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\')
            out << '\\' << char(c);
        else if (c < 32)
            out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c);
        else
            out << char(c);
    }
    return out.str() + '"';
}
Vec3 edge(const aeroacoustics::Node &n, bool leading) {
    Vec3 p = n.aero_center;
    for (int i = 0; i < 3; ++i)
        p[i] += n.section.chord * (-n.airfoil_reference[1] * n.global_to_local[i] +
                                   ((leading ? 0. : 1.) - n.airfoil_reference[0]) * n.global_to_local[3 + i]);
    return p;
}
} // namespace
struct EngineeringMetrics::Impl {
    MetricsOptions options;
    aeroacoustics::Parameters parameters;
    aeroacoustics::FrequencyBands bands;
    std::vector<Vec3> observers;
    Vec3 hub;
    std::optional<aeroacoustics::PropagationOptions> propagation;
    std::size_t nodes, values = 0;
    bool failed = false;
    std::vector<std::vector<ArrivalSeries>> paths, tones;
    std::vector<ArrivalSeries> wind;
    Impl(MetricsOptions o, aeroacoustics::Parameters p, std::vector<Vec3> receivers, std::size_t count,
         Vec3 h, std::optional<aeroacoustics::PropagationOptions> prop)
        : options(std::move(o)), parameters(std::move(p)),
          bands(aeroacoustics::FrequencyBands::openfast_reference(parameters.freqlist)),
          observers(std::move(receivers)), hub(h), propagation(std::move(prop)), nodes(count),
          paths(observers.size(), std::vector<ArrivalSeries>(2 * count)),
          tones(observers.size(), std::vector<ArrivalSeries>(options.tones.size())), wind(observers.size()) {
        options.validate();
        if (!nodes || observers.empty())
            throw std::invalid_argument("Metrics need sources and receivers");
        if (options.apparent_power && propagation)
            throw std::invalid_argument(
                "ApparentPower requires the reference free field; cannot invert propagated receiver levels");
        if (propagation && options.retarded_time && propagation->ground != aeroacoustics::GroundModel::none)
            throw std::invalid_argument(
                "RetardedTime currently requires direct paths: coherent ground image has a different delay");
        if (propagation && options.retarded_time &&
            std::abs(propagation->sound_speed - parameters.spdsound) > 1e-8)
            throw std::invalid_argument("RetardedTime requires equal source and propagation sound speeds");
        if (propagation && !options.tones.empty() &&
            (propagation->ground != aeroacoustics::GroundModel::none || !propagation->screens.empty()))
            throw std::invalid_argument("Tones currently support free field and air absorption only");
    }
};
EngineeringMetrics::EngineeringMetrics(MetricsOptions o, aeroacoustics::Parameters p, std::vector<Vec3> obs,
                                       std::size_t n, Vec3 h,
                                       std::optional<aeroacoustics::PropagationOptions> prop)
    : impl_(std::make_unique<Impl>(std::move(o), std::move(p), std::move(obs), n, h, std::move(prop))) {}
EngineeringMetrics::~EngineeringMetrics() = default;
EngineeringMetrics::EngineeringMetrics(const EngineeringMetrics &o)
    : impl_(std::make_unique<Impl>(*o.impl_)) {}
EngineeringMetrics &EngineeringMetrics::operator=(const EngineeringMetrics &o) {
    if (this != &o)
        impl_ = std::make_unique<Impl>(*o.impl_);
    return *this;
}
void EngineeringMetrics::append(double time, double hub_wind, double rotor_speed,
                                const std::vector<std::vector<aeroacoustics::Node>> &blades,
                                std::size_t first_node, std::size_t first_observer,
                                const aeroacoustics::Snapshot &block) {
    auto &s = *impl_;
    if (s.failed)
        throw std::logic_error("Metrics history failed; restore a valid copy or reconstruct");
    s.failed = true;
    const auto nf = s.bands.size();
    if (!std::isfinite(hub_wind) || hub_wind < 0 || hub_wind / s.options.wind_bin_width > 1000000)
        throw std::invalid_argument("Hub wind exceeds supported bin range");
    if (first_observer + block.size() > s.observers.size() || blades.empty() ||
        first_node >= blades.front().size())
        throw std::invalid_argument("Metrics source/observer shape");
    const auto per_blade = blades.front().size() - first_node;
    if (per_blade * blades.size() != s.nodes)
        throw std::invalid_argument("Metrics selected node count");
    const std::size_t increment = block.size() * (s.nodes * 2 * (nf + 5) + s.options.tones.size() * 6 + 6);
    if (increment > s.options.max_values - s.values)
        throw std::length_error(
            "Metrics history exceeds MaxValues; shorten duration/map or raise the explicit limit");
    s.values += increment;
    const double c = s.parameters.spdsound;
    for (std::size_t local = 0; local < block.size(); ++local) {
        const auto o = first_observer + local;
        const auto &receiver = s.observers[o];
        if (block[local].size() != s.nodes)
            throw std::invalid_argument("Metrics snapshot shape");
        for (std::size_t n = 0; n < s.nodes; ++n) {
            const auto &node = blades.at(n / per_blade).at(first_node + n % per_blade);
            for (int side = 0; side < 2; ++side) {
                std::vector<double> power(nf, 0.);
                for (std::size_t m = 0; m < aeroacoustics::mechanism_count; ++m) {
                    if ((m == aeroacoustics::index(aeroacoustics::Mechanism::inflow)) != bool(side))
                        continue;
                    if (block[local][n][m].size() != nf)
                        throw std::invalid_argument("Metrics frequency shape");
                    for (std::size_t f = 0; f < nf; ++f)
                        power[f] += aeroacoustics::relative_power(block[local][n][m][f]);
                }
                double arrival =
                    time + (s.options.retarded_time ? distance(edge(node, side == 1), receiver) / c : 0.);
                s.paths[o][2 * n + side].append({arrival, time, hub_wind, std::move(power)});
            }
        }
        s.wind[o].append(
            {time + (s.options.retarded_time ? distance(s.hub, receiver) / c : 0.), time, hub_wind, {0.}});
        for (std::size_t k = 0; k < s.options.tones.size(); ++k) {
            const auto &tone = s.options.tones[k];
            const auto &node = blades.at(tone.blade).at(tone.node);
            double r = distance(node.aero_center, receiver),
                   frequency = tone.frequency_hz + tone.rotor_order * rotor_speed / (2 * pi);
            if (r <= 0 || frequency <= 0)
                throw std::invalid_argument("Invalid tone distance/frequency");
            double db = tone.level_db + 20 * std::log10(tone.reference_distance / r);
            s.tones[o][k].append({time + (s.options.retarded_time ? r / c : 0.),
                                  time,
                                  hub_wind,
                                  {aeroacoustics::relative_power(db)},
                                  frequency,
                                  r});
        }
    }
    s.failed = false;
}
void EngineeringMetrics::write(const std::filesystem::path &directory) const {
    const auto &s = *impl_;
    if (s.failed)
        throw std::logic_error("Cannot write failed metrics history");
    const auto &opt = s.options;
    const auto nf = s.bands.size();
    double begin = opt.start, end = opt.end, source_dt = 0;
    auto intersect = [&](const ArrivalSeries &path) {
        const auto &a = path.samples();
        if (a.size() < 2)
            throw std::runtime_error("Metrics require at least two acoustic samples");
        begin = std::max(begin, a.front().reception);
        end = std::min(end, a.back().reception);
        source_dt = std::max(source_dt, a[1].emission - a[0].emission);
    };
    for (std::size_t o = 0; o < s.observers.size(); ++o) {
        for (const auto &p : s.paths[o])
            intersect(p);
        for (const auto &p : s.tones[o])
            intersect(p);
        intersect(s.wind[o]);
    }
    if (end - begin < opt.receiver_dt)
        throw std::runtime_error(
            "No common receiver-time window; extend source history or change metrics interval");
    if (source_dt * opt.am_max_hz >= .5)
        throw std::invalid_argument("Source acoustic sampling cannot resolve requested AM frequencies");
    std::vector<double> time;
    const auto count = std::size_t(std::floor((end - begin) / opt.receiver_dt + 1e-9));
    for (std::size_t i = 0; i <= count; ++i)
        time.push_back(std::min(end, begin + i * opt.receiver_dt));
    if (end - time.back() > 1e-10)
        time.push_back(end);
    diagnostics::CheckedOutput history, map, bins, am, tonal, meta;
    history.open(directory / "receiver_history.csv");
    map.open(directory / "receiver_map.csv");
    bins.open(directory / "wind_bins.csv");
    am.open(directory / "am_windows.csv");
    tonal.open(directory / "receiver_tones.csv");
    meta.open(directory / "metrics.json");
    history << std::setprecision(17)
            << "observer,time_receiver_s,hub_wind_m_s,turbine_LA_dB,total_with_background_LA_dB";
    for (double f : s.parameters.freqlist)
        history << ",band_" << f << "_dB";
    history << '\n';
    map << std::setprecision(17)
        << "observer,x_m,y_m,z_m,start_receiver_s,end_receiver_s,duration_s,LAeq_turbine_dB,LAeq_with_"
           "background_dB,L5_dBA,L50_dBA,L95_dBA,apparent_LWA_freefield_dB\n";
    bins << std::setprecision(17) << "observer,wind_lower_m_s,wind_upper_m_s,duration_s,LAeq_turbine_dB\n";
    am << std::setprecision(17)
       << "observer,start_receiver_s,end_receiver_s,resolved,peak_modulation_Hz,harmonic_peak_to_peak_dB,"
          "reconstructed_P95_minus_P5_dB\n";
    tonal << std::setprecision(17)
          << "observer,tone,time_receiver_s,frequency_received_Hz,SPL_unweighted_dB,tone_to_broadband_band_"
             "dB\n";
    const auto weights = aeroacoustics::a_weighting(s.parameters.freqlist);
    const double background =
        opt.background_laeq_db ? aeroacoustics::relative_power(*opt.background_laeq_db) : 0.;
    for (std::size_t o = 0; o < s.observers.size(); ++o) {
        std::vector<double> energies, wind, levels;
        for (double t : time) {
            std::vector<double> bands(nf, 0.);
            for (const auto &p : s.paths[o]) {
                auto sample = p.at(t);
                for (std::size_t f = 0; f < nf; ++f)
                    bands[f] += sample.power[f];
            }
            double total = 0;
            for (std::size_t f = 0; f < nf; ++f)
                total += bands[f] * (s.parameters.aweighting ? 1 : aeroacoustics::relative_power(weights[f]));
            const auto broadband = bands;
            for (std::size_t k = 0; k < s.tones[o].size(); ++k) {
                const auto sample = s.tones[o][k].at(t);
                const double frequency = sample.frequency_hz * s.tones[o][k].doppler(t);
                double p = sample.power[0];
                if (s.propagation && s.propagation->absorption)
                    p *= aeroacoustics::relative_power(
                        -s.propagation->atmosphere.absorption_db_per_m(frequency) * sample.distance_m);
                std::size_t band = 0;
                while (band < nf && !(frequency >= s.bands.values()[band].lower_hz &&
                                      frequency < s.bands.values()[band].upper_hz))
                    ++band;
                if (band == nf)
                    throw std::out_of_range("Received tone falls outside selected acoustic bands");
                const double gain = aeroacoustics::relative_power(aeroacoustics::a_weighting({frequency})[0]);
                const double broad =
                    broadband[band] /
                    (s.parameters.aweighting ? aeroacoustics::relative_power(weights[band]) : 1.);
                tonal << o + 1 << ',' << k + 1 << ',' << t << ',' << frequency << ',' << level(p) << ',';
                if (broad > 0)
                    tonal << level(p) - level(broad); // empty: no broadband denominator
                tonal << '\n';
                total += p * gain;
                bands[band] += p * (s.parameters.aweighting ? gain : 1.);
            }
            const double w = s.wind[o].at(t).wind;
            wind.push_back(w);
            energies.push_back(total);
            levels.push_back(level(total));
            history << o + 1 << ',' << t << ',' << w << ',' << level(total) << ','
                    << level(total + background);
            for (double p : bands)
                history << ',' << level(p);
            history << '\n';
        }
        auto stats = level_statistics(time, energies);
        auto with_background = energies;
        for (auto &p : with_background)
            p += background;
        const auto &position = s.observers[o];
        map << o + 1 << ',' << position[0] << ',' << position[1] << ',' << position[2] << ',' << begin << ','
            << end << ',' << stats.duration << ',' << stats.leq_db << ','
            << level_statistics(time, with_background).leq_db << ',' << stats.l5_db << ',' << stats.l50_db
            << ',' << stats.l95_db << ',';
        if (opt.apparent_power)
            map << apparent_sound_power(stats.leq_db, distance(s.hub, position));
        map << '\n';
        // Split every time interval at wind-bin boundaries; integrate the linear energy exactly.
        std::map<long long, std::pair<double, double>> accumulated;
        for (std::size_t i = 1; i < time.size(); ++i) {
            std::vector<double> cuts{0, 1};
            double w0 = wind[i - 1], w1 = wind[i], delta = w1 - w0;
            if (delta != 0) {
                long long first =
                              static_cast<long long>(std::floor(std::min(w0, w1) / opt.wind_bin_width)) + 1,
                          last = static_cast<long long>(std::floor(std::max(w0, w1) / opt.wind_bin_width));
                if (last - first > 100000)
                    throw std::length_error("Excessive wind bins");
                for (auto k = first; k <= last; ++k) {
                    double f = (k * opt.wind_bin_width - w0) / delta;
                    if (f > 0 && f < 1)
                        cuts.push_back(f);
                }
            }
            std::sort(cuts.begin(), cuts.end());
            for (std::size_t j = 1; j < cuts.size(); ++j) {
                double a = cuts[j - 1], b = cuts[j], dt = (b - a) * (time[i] - time[i - 1]);
                auto bin =
                    static_cast<long long>(std::floor((w0 + .5 * (a + b) * delta) / opt.wind_bin_width));
                auto &v = accumulated[bin];
                v.first += dt;
                v.second += dt * (energies[i - 1] + .5 * (a + b) * (energies[i] - energies[i - 1]));
            }
        }
        for (const auto &[bin, value] : accumulated)
            bins << o + 1 << ',' << bin * opt.wind_bin_width << ',' << (bin + 1) * opt.wind_bin_width << ','
                 << value.first << ',' << level(value.second / value.first) << '\n';
        const auto n = std::size_t(std::llround(opt.am_window / opt.receiver_dt));
        for (std::size_t first = 0; first + n < time.size(); first += n) {
            std::vector<double> segment(levels.begin() + first, levels.begin() + first + n);
            auto result = modulation(segment, opt.receiver_dt, opt.am_min_hz, opt.am_max_hz);
            am << o + 1 << ',' << time[first] << ',' << time[first] + n * opt.receiver_dt << ','
               << result.resolved << ',' << result.frequency_hz << ',' << result.harmonic_depth_db << ','
               << result.percentile_depth_db << '\n';
        }
    }
    history.finish();
    map.finish();
    bins.finish();
    am.finish();
    tonal.finish();
    meta << std::setprecision(17) << "{\n  \"requested_start_s\":" << opt.start
         << ",\n  \"requested_end_s\":" << opt.end << ",\n  \"actual_start_receiver_s\":" << begin
         << ",\n  \"actual_end_receiver_s\":" << end << ",\n  \"receiver_dt_s\":" << opt.receiver_dt
         << ",\n  \"source_dt_s\":" << source_dt
         << ",\n  \"retarded_time\":" << (opt.retarded_time ? "true" : "false")
         << ",\n  \"band_weighting\":" << quote(s.parameters.aweighting ? "A" : "unweighted")
         << ",\n  \"statistics_weighting\":\"A\",\n  \"reference_pressure_Pa\":0.00002,\n  "
            "\"time_average\":\"piecewise-linear mean-square pressure; exact trapezoidal integral\",\n  "
            "\"percentiles\":\"duration-weighted exceedance levels of interpolated energy\",\n"
         << "  \"wind_binning\":\"horizontal hub speed at source time delayed from initial hub position; "
            "exact interval splitting\",\n  \"background\":\"independent constant A-weighted background "
            "added to turbine prediction; no measured background subtraction\",\n  \"background_LAeq_dB\":";
    if (opt.background_laeq_db)
        meta << *opt.background_laeq_db;
    else
        meta << "null";
    meta << ",\n  \"am_method\":\"descriptive three-harmonic DFT reconstruction; not IOA reference or IEC "
            "compliance\",\n  \"am_window_s\":"
         << opt.am_window << ",\n  \"am_frequency_range_Hz\":[" << opt.am_min_hz << ',' << opt.am_max_hz
         << "],\n  \"partial_am_windows\":\"excluded\",\n  \"broadband_motion\":\"separate leading/trailing "
            "edge delay of level envelopes; no extra Doppler frequency remap or amplitude multiplier\",\n"
         << "  \"tones\":\"independent incoherent isotropic input lines; Doppler from arrival-time "
            "derivative; no audio phase or mechanical source prediction\",\n  \"tone_metric\":\"line energy "
            "relative to containing broadband band; not IEC tonal audibility\",\n"
         << "  \"apparent_power\":"
         << quote(opt.apparent_power ? "directional equivalent 4*pi*r^2 using initial hub distance; "
                                       "reference free field only, not IEC emission measurement"
                                     : "disabled")
         << ",\n  \"history_values\":" << s.values << ",\n  \"tone_sources\":[";
    for (std::size_t k = 0; k < opt.tones.size(); ++k) {
        if (k)
            meta << ',';
        const auto &t = opt.tones[k];
        meta << "{\"name\":" << quote(t.name) << ",\"provenance\":" << quote(t.provenance)
             << ",\"blade\":" << t.blade + 1 << ",\"node\":" << t.node + 1
             << ",\"frequency_Hz\":" << t.frequency_hz << ",\"rotor_order\":" << t.rotor_order
             << ",\"reference_SPL_dB\":" << t.level_db << ",\"reference_distance_m\":" << t.reference_distance
             << '}';
    }
    meta << "]\n}\n";
    meta.finish();
}
} // namespace turbine
