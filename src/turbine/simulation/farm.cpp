#include "turbine/farm.hpp"
#include "acoustic_levels.hpp"
#include "turbine/file_output.hpp"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <set>
namespace turbine {
namespace {
void require(bool ok, const char *s) {
    if (!ok)
        throw std::invalid_argument(s);
}
bool near(double a, double b) { return std::abs(a - b) < 1e-9 * std::max({1., std::abs(a), std::abs(b)}); }
std::string csv(const std::string &s) {
    std::string r = "\"";
    for (char c : s) {
        if (c == '"')
            r += '"';
        if (c == '\n' || c == '\r')
            r += ' ';
        else
            r += c;
    }
    return r + '"';
}
std::string json(const std::string &s) {
    std::string r = "\"";
    for (unsigned char c : s) {
        if (c == '"' || c == '\\')
            r += '\\';
        if (c < 32)
            r += ' ';
        else
            r += char(c);
    }
    return r + '"';
}
double db(double p) {
    require(std::isfinite(p) && p >= 0, "Non-finite/negative farm power");
    return p ? 10 * std::log10(p) : -INFINITY;
}
LevelStatistics stats(const std::vector<double> &time, const std::vector<double> &power, double start,
                      double end) {
    require(time.size() == power.size() && time.size() >= 2 && start >= time.front() - 1e-9 &&
                end <= time.back() + 1e-9,
            "Farm statistics window outside acoustic samples");
    auto at = [&](double t) {
        auto it = std::upper_bound(time.begin(), time.end(), t);
        const auto j = std::clamp<std::size_t>(it - time.begin(), 1, time.size() - 1);
        const double q = std::clamp((t - time[j - 1]) / (time[j] - time[j - 1]), 0., 1.);
        return power[j - 1] + q * (power[j] - power[j - 1]);
    };
    std::vector<double> t{start}, p{at(start)};
    for (std::size_t i = 0; i < time.size(); ++i)
        if (time[i] > start && time[i] < end) {
            t.push_back(time[i]);
            p.push_back(power[i]);
        }
    t.push_back(end);
    p.push_back(at(end));
    return level_statistics(t, p);
}
void alignment(const BladeStructure &s, const Vec3 &direction) {
    const Vec3 horizontal{s.shaft[0], s.shaft[1], 0};
    require(norm(unit(horizontal) - direction) < 1e-7 && dot(s.shaft, direction) > std::cos(10 * deg),
            "Jensen farm requires aligned horizontal yaw and tilt below 10 degrees");
}
} // namespace
void run_farm(const FarmOptions &o, const std::filesystem::path &output) {
    require(!o.provenance.empty() && !o.turbines.empty() && o.turbines.size() <= 1000 && !o.observers.empty(),
            "Missing farm metadata/turbines/receivers");
    require(std::isfinite(o.duration) && std::isfinite(o.statistics_start) &&
                o.duration > o.statistics_start && o.statistics_start >= 0,
            "Invalid farm duration/statistics window");
    require(std::isfinite(o.expansion) && o.expansion > 0 && std::isfinite(o.max_ct) && o.max_ct > 0 &&
                o.max_ct < 1 && o.max_values > 0,
            "Invalid wake expansion/CT/memory limit");
    for (const auto &p : o.observers)
        for (double x : p)
            require(std::isfinite(x), "Invalid farm receiver");
    if (o.propagation) {
        o.propagation->validate();
        require(o.propagation->ground == aeroacoustics::GroundModel::none && o.propagation->screens.empty(),
                "Farm currently supports free field/air absorption only");
    }
    std::vector<Case> cases;
    std::vector<Vec3> hubs;
    std::vector<double> radii;
    std::set<std::string> names;
    for (const auto &u : o.turbines) {
        require(!u.name.empty() && names.insert(u.name).second, "Duplicate/empty turbine name");
        for (double v : u.origin)
            require(std::isfinite(v), "Invalid turbine origin");
        require(u.origin[2] == 0, "Farm requires a flat site");
        if (u.tower)
            u.tower->validate();
        cases.emplace_back(u.case_file);
        BladeStructure structure(cases.back());
        hubs.push_back(u.origin + structure.hub);
        radii.push_back(structure.length + structure.hub_radius);
    }
    const auto &reference = cases.front();
    const AcousticConfiguration config(reference);
    const auto bands = aeroacoustics::FrequencyBands::openfast_reference(config.parameters.freqlist);
    const std::size_t nf = bands.size(), no = o.observers.size();
    const auto &ambient = reference.wind;
    require(ambient.speed > 0 && ambient.upflow == 0, "Farm requires positive horizontal steady wind");
    const Vec3 direction{std::cos(ambient.propagation), -std::sin(ambient.propagation), 0};
    for (const auto &c : cases) {
        AcousticConfiguration a(c);
        alignment(BladeStructure(c), direction);
        require(near(c.wind.speed, ambient.speed) &&
                    near(c.wind.reference_height, ambient.reference_height) &&
                    near(c.wind.exponent, ambient.exponent) &&
                    near(c.wind.propagation, ambient.propagation) && c.wind.upflow == 0,
                "Farm turbines must share ambient wind definition");
        require(near(c.rho, reference.rho) && near(c.sound_speed, reference.sound_speed),
                "Farm medium mismatch");
        require(near(c.dt, reference.dt) && near(a.sample_interval, config.sample_interval) &&
                    near(a.start, config.start) && a.parameters.freqlist == config.parameters.freqlist,
                "Farm time grid/frequency mismatch");
        require(c.acoustic.integer("NrOutFile") >= 2, "Farm requires NrOutFile >= 2 for band energies");
    }
    require(config.sample_interval > 0 &&
                near(o.duration / config.sample_interval, std::round(o.duration / config.sample_interval)) &&
                near(o.duration / reference.dt, std::round(o.duration / reference.dt)),
            "Farm duration must be a multiple of structural/acoustic steps");
    require(o.statistics_start >= config.start, "Farm statistics start precedes acoustic start");
    const long double estimate = (std::ceil(o.duration / reference.dt) + 1) * 3 * cases.size() +
                                 (std::ceil(o.duration / config.sample_interval) + 1) * (no * nf + no + 1) +
                                 no * nf * o.sources.size();
    require(estimate <= o.max_values, "Farm history exceeds MaxValues");
    std::vector<std::size_t> order(cases.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) {
        const double aa = dot(hubs[a], direction), bb = dot(hubs[b], direction);
        return aa == bb ? o.turbines[a].name < o.turbines[b].name : aa < bb;
    });
    for (std::size_t i = 0; i < order.size(); ++i)
        for (std::size_t j = i + 1; j < order.size(); ++j) {
            const auto a = order[i], b = order[j];
            const Vec3 delta = hubs[b] - hubs[a];
            const double x = dot(delta, direction), crosswind = norm(delta - x * direction);
            require(norm(delta) > radii[a] + radii[b], "Overlapping turbine exclusion spheres");
            if (o.wakes && crosswind <= radii[a] + o.expansion * x + radii[b])
                require(x >= 4 * radii[a],
                        "Interacting turbines must be separated by at least two upstream diameters");
        }
    std::vector<std::vector<std::vector<double>>> external;
    std::set<std::string> source_names;
    for (const auto &s : o.sources) {
        require(source_names.insert(s.name).second, "Duplicate stationary source name");
        std::vector<std::vector<double>> powers;
        for (const auto &r : o.observers)
            powers.push_back(s.receiver_power(r, bands, reference.rho, reference.sound_speed, o.propagation));
        external.push_back(std::move(powers));
    }
    require(std::filesystem::create_directory(output),
            "Farm output directory must be new; parent must exist");
    diagnostics::CheckedOutput thrust, contributions, manifest;
    thrust.open(output / "wake_diagnostics.csv");
    contributions.open(output / "source_contributions.csv");
    manifest.open(output / "stationary_sources.csv");
    thrust << std::setprecision(17) << "turbine,time_s,thrust_N,hub_inflow_m_s,ct_raw,ct_used,ct_limited\n";
    contributions << std::setprecision(17) << "category,source,observer,start_s,end_s,LAeq_dBA\n";
    manifest << std::setprecision(17)
             << "name,kind,x_m,y_m,z_m,axis_x,axis_y,axis_z,directivity,frequency_hz,Lw_unweighted_dB,"
                "provenance\n";
    for (const auto &s : o.sources)
        for (std::size_t f = 0; f < s.frequencies.size(); ++f)
            manifest << csv(s.name) << ',' << s.kind << ',' << s.position[0] << ',' << s.position[1] << ','
                     << s.position[2] << ',' << s.axis[0] << ',' << s.axis[1] << ',' << s.axis[2] << ','
                     << s.directivity << ',' << s.frequencies[f] << ',' << s.levels_db[f] << ','
                     << csv(s.provenance) << '\n';
    std::vector<std::shared_ptr<const WakeHistory>> histories;
    std::vector<double> times;
    std::vector<std::vector<double>> total;
    std::size_t limited = 0;
    for (std::size_t sequence = 0; sequence < order.size(); ++sequence) {
        const auto index = order[sequence];
        const auto &u = o.turbines[index];
        auto field = std::make_shared<const FarmWind>(
            ambient, u.origin, o.expansion,
            o.wakes ? histories : std::vector<std::shared_ptr<const WakeHistory>>{});
        RunOptions options;
        options.duration = o.duration;
        options.solver.wind = field;
        options.solver.tower = u.tower;
        options.solver.controller = u.controller;
        options.surfaces = u.surfaces;
        options.propagation = o.propagation;
        for (const auto &r : o.observers)
            options.observers.push_back(r - u.origin);
        Simulation simulation(cases[index], options);
        FileOutput files(output / ("turbine_" + std::to_string(index + 1)));
        files.begin(simulation.layout());
        auto history = std::make_shared<WakeHistory>();
        history->hub = hubs[index];
        history->radius = radii[index];
        std::size_t sample = 0;
        std::vector<std::vector<double>> individual(no);
        const auto weights = aeroacoustics::a_weighting(config.parameters.freqlist);
        while (auto frame = simulation.next()) {
            files.write(*frame);
            const auto &rotor = simulation.solver().rotor();
            alignment(rotor.structure(), direction);
            require(norm(rotor.structure().hub + u.origin - hubs[index]) < 1e-7,
                    "Farm wake centers must remain fixed");
            const double speed = dot(field->at(frame->time, rotor.structure().hub), direction),
                         force = rotor.aerodynamic_thrust(simulation.solver().aerodynamic(), direction);
            const double raw =
                force / (.5 * reference.rho * pi * radii[index] * radii[index] * speed * speed);
            require(std::isfinite(raw), "Non-finite farm thrust coefficient");
            const double ct = std::clamp(raw, 0., o.max_ct);
            if (raw != ct)
                ++limited;
            history->samples.push_back({frame->time, speed, ct});
            thrust << csv(u.name) << ',' << frame->time << ',' << force << ',' << speed << ',' << raw << ','
                   << ct << ',' << (raw != ct ? 1 : 0) << '\n';
            if (frame->acoustics) {
                const auto &powers = frame->acoustics->power[1];
                require(powers.size() == no * nf, "Farm spectrum shape mismatch");
                if (sequence == 0) {
                    times.push_back(frame->time);
                    total.emplace_back(no * nf, 0);
                } else
                    require(sample < times.size() && near(times[sample], frame->time),
                            "Farm acoustic clocks differ");
                for (std::size_t r = 0; r < no; ++r) {
                    double sum = 0;
                    for (std::size_t f = 0; f < nf; ++f) {
                        const double p =
                            powers[r * nf + f] * (simulation.layout().parameters.aweighting
                                                      ? 1
                                                      : aeroacoustics::relative_power(weights[f]));
                        require(std::isfinite(p) && p >= 0, "Invalid turbine energy");
                        sum += p;
                        total[sample][r * nf + f] += p;
                    }
                    individual[r].push_back(sum);
                }
                ++sample;
            }
        }
        require(sample == times.size(), "Missing farm acoustic samples");
        files.finish(simulation.summary());
        history->validate();
        histories.push_back(history);
        for (std::size_t r = 0; r < no; ++r) {
            const auto s = stats(times, individual[r], o.statistics_start, o.duration);
            contributions << "blades," << csv(u.name) << ',' << r + 1 << ',' << o.statistics_start << ','
                          << o.duration << ',' << s.leq_db << '\n';
        }
    }
    for (std::size_t i = 0; i < external.size(); ++i)
        for (std::size_t r = 0; r < no; ++r) {
            const auto &p = external[i][r];
            const double sum = std::accumulate(p.begin(), p.end(), 0.);
            contributions << "stationary," << csv(o.sources[i].name) << ',' << r + 1 << ','
                          << o.statistics_start << ',' << o.duration << ',' << db(sum) << '\n';
            for (auto &sample : total)
                for (std::size_t f = 0; f < nf; ++f)
                    sample[r * nf + f] += p[f];
        }
    diagnostics::CheckedOutput acoustic, map;
    acoustic.open(output / "farm_history.csv");
    map.open(output / "farm_receivers.csv");
    acoustic << std::setprecision(17) << "time_source_s,observer,LA_dBA";
    for (double f : config.parameters.freqlist)
        acoustic << ",A_band_" << f << "_Hz_dB";
    acoustic << '\n';
    map << std::setprecision(17)
        << "observer,x_m,y_m,z_m,start_source_s,end_source_s,LAeq_dBA,L5_dBA,L50_dBA,L95_dBA\n";
    for (std::size_t r = 0; r < no; ++r) {
        std::vector<double> powers;
        for (std::size_t t = 0; t < times.size(); ++t) {
            double sum = 0;
            for (std::size_t f = 0; f < nf; ++f)
                sum += total[t][r * nf + f];
            powers.push_back(sum);
            acoustic << times[t] << ',' << r + 1 << ',' << db(sum);
            for (std::size_t f = 0; f < nf; ++f)
                acoustic << ',' << db(total[t][r * nf + f]);
            acoustic << '\n';
        }
        const auto s = stats(times, powers, o.statistics_start, o.duration);
        const auto &p = o.observers[r];
        map << r + 1 << ',' << p[0] << ',' << p[1] << ',' << p[2] << ',' << o.statistics_start << ','
            << o.duration << ',' << s.leq_db << ',' << s.l5_db << ',' << s.l50_db << ',' << s.l95_db << '\n';
    }
    acoustic.finish();
    map.finish();
    thrust.finish();
    contributions.finish();
    manifest.finish();
    diagnostics::CheckedOutput metadata;
    metadata.open(output / "farm.json");
    metadata << std::setprecision(17)
             << "{\n  \"solver\":\"independent C++ turbine ensemble\",\n  \"provenance\":"
             << json(o.provenance) << ",\n  \"wake_model\":"
             << json(o.wakes ? "quasi-steady Jensen top-hat; dimensional RSS deficits" : "disabled")
             << ",\n  \"wake_expansion\":" << o.expansion << ",\n  \"max_ct\":" << o.max_ct
             << ",\n  \"ct_limited_samples\":" << limited << ",\n  \"duration_s\":" << o.duration
             << ",\n  \"statistics_start_s\":" << o.statistics_start
             << ",\n  \"air_density_kg_m3\":" << reference.rho
             << ",\n  \"sound_speed_m_s\":" << reference.sound_speed
             << ",\n  \"air_absorption\":" << (o.propagation && o.propagation->absorption ? "true" : "false")
             << ",\n  \"time_basis\":\"common source time; no acoustic or wake transport delay\",\n  "
                "\"superposition\":\"incoherent mean-square pressure; no phase or cross terms\",\n  "
                "\"field_validated\":false,\n  \"turbines\":[";
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto j = order[i];
        const auto &u = o.turbines[j];
        if (i)
            metadata << ',';
        metadata << "{\"name\":" << json(u.name)
                 << ",\"directory\":" << json("turbine_" + std::to_string(j + 1))
                 << ",\"input\":" << json(std::filesystem::absolute(u.case_file).generic_string())
                 << ",\"origin\":[" << u.origin[0] << ',' << u.origin[1] << ",0]}";
    }
    metadata << "]\n}\n";
    metadata.finish();
}
} // namespace turbine
