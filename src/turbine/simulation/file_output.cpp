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
std::string json_text(const std::string &s) {
    constexpr char hex[] = "0123456789abcdef";
    std::string result = "\"";
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') {
            result += '\\';
            result += c;
        } else if (c < 32) {
            result += "\\u00";
            result += hex[c >> 4];
            result += hex[c & 15];
        } else
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
    if (layout.output_count < 1 || layout.output_count > 4 || !layout.blades)
        throw std::invalid_argument("Invalid output layout");
    if (!layout.acoustic_metadata)
        throw std::invalid_argument("Missing acoustic metadata");
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
            outputs_[k] << "\t(" << layout.acoustic_metadata->level_unit() << ')';
        outputs_[k] << '\n' << std::setprecision(12);
    }
    dynamics_.open(directory_ / "dynamics.csv");
    lookup_output_.open(directory_ / "lookup_diagnostics.csv");
    dynamics_ << std::setprecision(17) << "time";
    for (const auto &channel : layout.dynamics)
        dynamics_ << ',' << channel.name;
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
    if (frame.generalized) {
        if (frame.generalized->q.size() != layout_->dofs.size() ||
            frame.generalized->qd.size() != layout_->dofs.size())
            throw std::invalid_argument("Generalized output state shape mismatch");
        for (const auto &channel : layout_->dynamics) {
            const auto &values = channel.velocity ? frame.generalized->qd : frame.generalized->q;
            dynamics_ << ',' << values.at(channel.dof);
        }
    } else {
        // Compatibility for callers constructing the original three-field StepView.
        if (layout_->dofs.size() != FixedBaseBladeBackend::blades * FixedBaseBladeBackend::modes_per_blade)
            throw std::invalid_argument("Legacy state view requires the fixed-base layout");
        for (const auto &channel : layout_->dynamics) {
            const auto &state = frame.state.at(channel.dof / FixedBaseBladeBackend::modes_per_blade);
            dynamics_ << ','
                      << (channel.velocity ? state.qd : state.q)
                             .at(channel.dof % FixedBaseBladeBackend::modes_per_blade);
        }
    }
    dynamics_ << '\n';

    if (frame.acoustics) {
        const auto &metadata = frame.acoustics->metadata;
        if (!metadata || metadata->weighting != layout_->acoustic_metadata->weighting ||
            !metadata->bands.same_as(layout_->acoustic_metadata->bands))
            throw std::invalid_argument("Output acoustic quantity/weighting/bands mismatch");
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
              << ",\n  \"lookup_report\": \"lookup_diagnostics.csv\",\n  \"acoustic_metadata\": {\n";
    const auto &acoustic = *layout_->acoustic_metadata;
    metadata_
        << "    \"linear_quantity\": \"" << acoustic.linear_quantity()
        << "\",\n    \"level_quantity\": \"sound_pressure_level\",\n    \"reference_pressure_Pa\": 0.00002,\n"
        << "    \"weighting\": \"" << (acoustic.weighting == aeroacoustics::Weighting::a ? "A" : "unweighted")
        << "\",\n    \"frequency_unit\": \"Hz\",\n"
        << "    \"band_definition\": \"binary third-octave edges; OpenFAST nominal evaluation centers\",\n"
        << "    \"tno_bandwidth_convention\": \"OpenFAST legacy multiplier; not SI Hz bandwidth\",\n"
        << "    \"bands_center_lower_upper_Hz\": [";
    bool comma = false;
    for (const auto &b : acoustic.bands.values()) {
        if (comma)
            metadata_ << ',';
        comma = true;
        metadata_ << '[' << b.center_hz << ',' << b.lower_hz << ',' << b.upper_hz << ']';
    }
    metadata_ << "]\n  },\n  \"module_configuration\": {\"profile\": " << json_text(layout_->module_profile)
              << ", \"dof_count\": " << layout_->dofs.size()
              << ", \"integrator\": \"generalized_alpha\", \"dofs\": [";
    comma = false;
    for (const auto &d : layout_->dofs) {
        if (comma)
            metadata_ << ',';
        comma = true;
        metadata_ << "{\"name\":" << json_text(d.name) << ",\"unit\":" << json_text(d.unit)
                  << ",\"acceleration_scale\":" << d.acceleration_scale << '}';
    }
    metadata_ << "], \"coupling_blocks\": [";
    comma = false;
    for (const auto &b : layout_->coupling_blocks) {
        if (comma)
            metadata_ << ',';
        comma = true;
        metadata_ << "{\"name\":" << json_text(b.name) << ",\"offset\":" << b.offset << ",\"size\":" << b.size
                  << '}';
    }
    metadata_ << "]}\n}\n";
    metadata_.finish();

    finished_ = true;
}
} // namespace turbine
