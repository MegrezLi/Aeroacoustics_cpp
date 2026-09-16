#include "turbine/file_output.hpp"
#include "acoustic_levels.hpp"
#include <iomanip>
namespace turbine {
namespace {
std::string csv_text(const std::string &s) {
    std::string result = "\"";
    for (char c : s) {
        result += c;
        if (c == '"')
            result += c;
    }
    return result + '"';
}

} // namespace
FileOutput::FileOutput(std::filesystem::path directory) : directory_(std::move(directory)) {
    std::filesystem::create_directories(directory_);
    metadata_.open(directory_ / "run.json");
}
void FileOutput::begin(const OutputLayout &layout) {
    if (failed_)
        throw std::logic_error("Output failed; create a new sink before continuing");
    try {
        begin_output(layout);
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void FileOutput::begin_output(const OutputLayout &layout) {
    if (layout_)
        throw std::logic_error("Output already initialized");
    if (layout.output_count < 1 || layout.output_count > 4 || layout.blades != 3)
        throw std::invalid_argument("Invalid output layout");
    layout_ = layout;
    const auto &labels = layout.labels;
    const auto &prefix = layout.prefix;
    const auto output_count = layout.output_count;
    for (int k = 0; k < output_count; ++k) {
        outputs_[k].open(directory_ / (prefix + std::to_string(k + 1) + ".out"));
        masks_[k].open(directory_ / (prefix + std::to_string(k + 1) + ".mask"));
        masks_[k] << "Energy mask: 1=positive energy; 0=zero energy (legacy 0 dB placeholder)\nTime";
        for (const auto &label : labels[k])
            masks_[k] << '\t' << label;
        masks_[k] << "\n(s)";
        for (std::size_t j = 0; j < labels[k].size(); ++j)
            masks_[k] << "\t(flag)";
        masks_[k] << '\n' << std::setprecision(12);
        outputs_[k] << "C++ turbine aeroacoustics; reference sound pressure 20 uPa\nTime";
        for (const auto &label : labels[k])
            outputs_[k] << '\t' << label;
        outputs_[k] << "\n(s)";
        for (std::size_t j = 0; j < labels[k].size(); ++j)
            outputs_[k] << "\t("
                        << (k == 2
                                ? aeroacoustics::mechanism_registry[j % aeroacoustics::mechanism_count].unit
                                : "dB")
                        << ')';
        outputs_[k] << '\n' << std::setprecision(12);
    }
    dynamics_.open(directory_ / "dynamics.csv");
    lookup_output_.open(directory_ / "lookup_diagnostics.csv");
    dynamics_ << std::setprecision(17) << "time";
    for (int b = 1; b <= 3; ++b)
        dynamics_ << ",q" << b << "_flap1,q" << b << "_flap2,q" << b << "_edge,qd" << b << "_flap1,qd" << b
                  << "_flap2,qd" << b << "_edge";
    dynamics_ << '\n';
}
void FileOutput::write(const StepView &frame) {
    if (failed_)
        throw std::logic_error("Output failed; create a new sink before continuing");
    try {
        write_step(frame);
        last_time_ = frame.time;
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void FileOutput::write_step(const StepView &frame) {
    if (!layout_ || finished_)
        throw std::logic_error("Output is not writable");
    if (!std::isfinite(frame.time) || frame.time <= last_time_)
        throw std::invalid_argument("Output times must be finite and increasing");
    const auto output_count = layout_->output_count;
    const auto &labels = layout_->labels;
    dynamics_ << frame.time;
    for (const auto &s : frame.state)
        for (const auto &v : {s.q, s.qd})
            for (double x : v)
                dynamics_ << ',' << x;
    dynamics_ << '\n';

    if (frame.acoustics) {
        const auto &values = frame.acoustics->power;
        for (int k = 0; k < output_count; ++k)
            if (values[k].size() != labels[k].size())
                throw std::invalid_argument("Output channel count mismatch");
        for (int k = 0; k < output_count; ++k) {
            outputs_[k] << frame.time;
            masks_[k] << frame.time;
            for (std::size_t channel = 0; channel < values[k].size(); ++channel) {
                const double p = values[k][channel];
                double db;
                try {
                    db = aeroacoustics::output_decibels(p);
                } catch (const std::exception &e) {
                    throw std::runtime_error(std::string(e.what()) + " time=" + std::to_string(frame.time) +
                                             " file=" + std::to_string(k + 1) +
                                             " channel=" + labels[k][channel]);
                }
                outputs_[k] << '\t' << db;
                masks_[k] << '\t' << (p > 0 ? 1 : 0);
            }
            outputs_[k] << '\n';
            masks_[k] << '\n';
        }
    }
}
void FileOutput::finish(const RunSummary &summary) {
    if (failed_)
        throw std::logic_error("Output failed; create a new sink before continuing");
    try {
        finish_output(summary);
    } catch (...) {
        failed_ = true;
        throw;
    }
}
void FileOutput::finish_output(const RunSummary &summary) {
    if (!layout_ || finished_)
        throw std::logic_error("Output is not active");
    const auto output_count = layout_->output_count;
    const auto &lookup = summary.lookup;
    lookup_output_
        << "table,axis,blade,node,stage,calls,first_time_s,last_time_s,min_value,max_value,lower,upper\n"
        << std::setprecision(17);
    for (const auto &entry : lookup.entries) {
        const auto &[table, axis, blade, node, stage] = entry.first;
        const auto &s = entry.second;
        lookup_output_ << csv_text(table) << ',' << axis << ',' << blade << ',' << node << ',' << stage << ','
                       << s.calls << ',' << s.first_time << ',' << s.last_time << ',' << s.minimum << ','
                       << s.maximum << ',' << s.lower << ',' << s.upper << '\n';
    }
    for (int k = 0; k < output_count; ++k) {
        outputs_[k].finish();
        masks_[k].finish();
    }
    dynamics_.finish();
    lookup_output_.finish();
    metadata_ << std::setprecision(17) << "{\n  \"solver\": \"standalone C++\",\n  \"dt\": " << summary.dt
              << ",\n  \"duration\": " << summary.duration << ",\n  \"steps\": " << summary.steps
              << ",\n  \"acoustic_samples\": " << summary.acoustic_samples
              << ",\n  \"elapsed_seconds\": " << summary.elapsed_seconds << ",\n  \"structural_mode\": \""
              << (summary.structure.mode == SolverMode::scaled ? "scaled" : "reference") << "\""
              << ",\n  \"structural_iterations\": " << summary.structure.iterations
              << ",\n  \"structural_max_iterations\": " << summary.structure.max_iterations
              << ",\n  \"structural_acceleration_evaluations\": "
              << summary.structure.acceleration_evaluations
              << ",\n  \"structural_jacobian_builds\": " << summary.structure.jacobian_builds
              << ",\n  \"structural_last_residual\": " << summary.structure.residual
              << ",\n  \"structural_last_scaled_residual\": " << summary.structure.scaled_residual
              << ",\n  \"zero_energy_encoding\": \"0 dB placeholder; see matching .mask (0=zero, "
                 "1=positive energy)\""
              << ",\n  \"lookup_policy\": \""
              << (lookup.policy == diagnostics::LookupPolicy::error ? "error" : "clamp")
              << "\",\n  \"lookup_out_of_range_calls\": " << lookup.calls
              << ",\n  \"lookup_report\": \"lookup_diagnostics.csv\"\n}\n";
    metadata_.finish();

    finished_ = true;
}
} // namespace turbine
