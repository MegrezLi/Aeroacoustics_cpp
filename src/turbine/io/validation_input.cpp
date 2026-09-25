#include "checked_output.hpp"
#include "turbine/validation.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace turbine::validation {
namespace {
using Row = std::map<std::string, std::string>;
using Rows = std::vector<Row>;
void require(bool ok, const std::string &s) {
    if (!ok)
        throw std::invalid_argument(s);
}
std::vector<std::string> fields(const std::string &s) {
    std::vector<std::string> result;
    std::string cell;
    bool quoted = false, closed = false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        require(static_cast<unsigned char>(c) >= 32 || c == '\t', "Control character in CSV");
        if (quoted) {
            if (c == '"') {
                if (i + 1 < s.size() && s[i + 1] == '"') {
                    cell += '"';
                    ++i;
                } else {
                    quoted = false;
                    closed = true;
                }
            } else
                cell += c;
        } else if (c == ',') {
            result.push_back(cell);
            cell.clear();
            closed = false;
        } else if (c == '"') {
            require(cell.empty() && !closed, "Malformed CSV quote");
            quoted = true;
        } else {
            require(!closed, "Characters after CSV closing quote");
            cell += c;
        }
    }
    require(!quoted, "Multiline/unterminated CSV field is unsupported");
    result.push_back(cell);
    return result;
}
Rows read(const std::filesystem::path &p, const std::vector<std::string> &schema, bool allow_empty = false,
          bool extra = false) {
    std::ifstream in(p);
    require(bool(in), "Cannot read " + p.string());
    std::string line;
    require(bool(std::getline(in, line)), "Missing CSV header: " + p.string());
    if (line.compare(0, 3, "\xef\xbb\xbf") == 0)
        line.erase(0, 3);
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    const auto names = fields(line);
    std::set<std::string> unique(names.begin(), names.end());
    require(unique.size() == names.size() && !unique.count(""), "Duplicate/empty CSV column");
    for (const auto &s : schema)
        require(unique.count(s), "Missing CSV column " + s);
    require(extra || names.size() == schema.size(), "Unexpected CSV column");
    Rows rows;
    std::size_t lineno = 1;
    while (std::getline(in, line)) {
        ++lineno;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        require(!line.empty() && line.size() <= 1048576 && rows.size() < 100000,
                "Empty/oversized CSV row or table");
        const auto cells = fields(line);
        require(cells.size() == names.size(), "CSV column count at line " + std::to_string(lineno));
        Row row;
        for (std::size_t i = 0; i < names.size(); ++i)
            row[names[i]] = cells[i];
        rows.push_back(std::move(row));
    }
    require(in.eof(), "CSV read failure");
    require(allow_empty || !rows.empty(), "Empty CSV data");
    return rows;
}
double number(const std::string &s) {
    std::size_t end = 0;
    double v = std::stod(s, &end);
    require(end == s.size() && std::isfinite(v), "Expected finite number: " + s);
    return v;
}
double num(const Row &r, const std::string &s) { return number(r.at(s)); }
std::string fmt(double v) {
    require(std::isfinite(v), "Non-finite report value");
    std::ostringstream s;
    s << std::setprecision(17) << v;
    return s.str();
}
std::string quote(const std::string &s) {
    std::string result = "\"";
    for (char c : s) {
        if (c == '"')
            result += '"';
        result += c;
    }
    return result + '"';
}
void write(const std::filesystem::path &p, const std::vector<std::string> &header,
           const std::vector<std::vector<std::string>> &rows) {
    diagnostics::CheckedOutput out;
    out.open(p);
    auto line = [&](const std::vector<std::string> &values) {
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i)
                out << ',';
            out << quote(values[i]);
        }
        out << '\n';
    };
    line(header);
    for (const auto &r : rows)
        line(r);
    out.finish();
}
void fresh(const std::filesystem::path &output) {
    require(std::filesystem::create_directory(output), "Output directory must be new (parent must exist)");
}
void finish(const std::filesystem::path &output, const std::string &method, double k, std::size_t external) {
    diagnostics::CheckedOutput out;
    out.open(output / "report.json");
    out << "{\n  \"schema\": 1,\n  \"method\": \"" << method
        << "\",\n  \"coverage_factor\": " << (k > 0 ? fmt(k) : "null")
        << ",\n  \"declared_external_validation_rows\": " << external
        << ",\n  \"field_accuracy_established\": false,\n  \"automatic_calibration\": false,\n"
           "  \"note\": \"Data provenance is caller-declared. No compliance or universal accuracy pass is "
           "inferred.\"\n}\n";
    out.finish();
}
const std::vector<std::string> schema = {
    "id",        "group",        "split",        "origin",        "provenance", "receiver", "quantity",
    "weighting", "frequency_hz", "bandwidth_hz", "signal",        "start_s",    "end_s",    "x_m",
    "y_m",       "z_m",          "wind_mps",     "direction_deg", "surface",    "level_db", "standard_u_db"};
const std::vector<std::string> numeric_context = {
    "frequency_hz", "bandwidth_hz", "start_s", "end_s", "x_m", "y_m", "z_m", "wind_mps", "direction_deg"};
void validate_record(const Row &r, bool simulation) {
    for (const auto &s : schema)
        require(s == "standard_u_db" || !r.at(s).empty(), "Empty record field " + s);
    require(r.at("split") == "calibration" || r.at("split") == "validation",
            "Split must be calibration or validation");
    const auto &origin = r.at("origin");
    require(simulation ? origin == "simulation"
                       : (origin == "synthetic" || origin == "field" || origin == "wind_tunnel"),
            "Invalid data origin");
    require(r.at("quantity") == "LAeq" || r.at("quantity") == "band_Leq", "Unsupported quantity");
    require(r.at("weighting") == "A" || r.at("weighting") == "Z", "Unsupported weighting");
    require(r.at("signal") == "turbine_only" || r.at("signal") == "total",
            "Signal must distinguish turbine-only/total");
    for (const auto &s : numeric_context)
        num(r, s);
    num(r, "level_db");
    if (!r.at("standard_u_db").empty())
        require(num(r, "standard_u_db") >= 0, "Negative standard uncertainty");
    require(num(r, "start_s") >= 0 && num(r, "end_s") > num(r, "start_s"), "Invalid averaging interval");
    require(num(r, "wind_mps") >= 0 && num(r, "wind_mps") < 1000 && num(r, "direction_deg") >= 0 &&
                num(r, "direction_deg") < 360,
            "Invalid wind/direction");
    require(r.at("quantity") == "LAeq"
                ? (r.at("weighting") == "A" && num(r, "frequency_hz") == 0 && num(r, "bandwidth_hz") == 0)
                : (num(r, "frequency_hz") > 0 && num(r, "bandwidth_hz") > 0 &&
                   num(r, "bandwidth_hz") < 2 * num(r, "frequency_hz")),
            "Quantity/frequency/weighting mismatch");
}
bool near(double a, double b) { return std::abs(a - b) <= 1e-8 * std::max({1., std::abs(a), std::abs(b)}); }
} // namespace
void compare_files(const std::filesystem::path &observed, const std::filesystem::path &predicted,
                   const std::filesystem::path &output, double k) {
    require(std::isfinite(k) && k > 0, "Invalid coverage factor");
    const auto obs = read(observed, schema), pred = read(predicted, schema);
    require(obs.size() == pred.size(), "Observation/prediction counts differ; no silent dropping");
    std::map<std::string, const Row *> lookup;
    std::map<std::string, std::string> splits;
    std::set<std::string> ids;
    for (const auto &r : pred) {
        validate_record(r, true);
        require(lookup.emplace(r.at("id"), &r).second, "Duplicate prediction ID");
    }
    std::map<std::vector<std::string>, std::vector<Residual>> groups;
    std::vector<std::vector<std::string>> pairs, summaries;
    std::size_t external = 0;
    for (const auto &r : obs) {
        validate_record(r, false);
        require(ids.insert(r.at("id")).second, "Duplicate observation ID");
        const auto split = splits.emplace(r.at("group"), r.at("split"));
        require(split.second || split.first->second == r.at("split"), "Calibration/validation group leakage");
        require(lookup.count(r.at("id")), "Missing prediction ID " + r.at("id"));
        const auto &p = *lookup.at(r.at("id"));
        for (const auto &name : {"group", "split", "receiver", "quantity", "weighting", "signal", "surface"})
            require(r.at(name) == p.at(name), std::string("Context mismatch: ") + name);
        for (const auto &name : numeric_context)
            require(near(num(r, name), num(p, name)), "Context mismatch: " + name);
        const bool uncertainty_known = !r.at("standard_u_db").empty() && !p.at("standard_u_db").empty();
        const auto e =
            residual(num(p, "level_db"), num(r, "level_db"), uncertainty_known ? num(p, "standard_u_db") : 0,
                     uncertainty_known ? num(r, "standard_u_db") : 0);
        pairs.push_back({r.at("id"), r.at("group"), r.at("split"), r.at("origin"), fmt(e.error_db),
                         uncertainty_known ? fmt(e.combined_u_db) : "",
                         e.normalized_available ? fmt(e.normalized_error) : "",
                         e.normalized_available ? (std::abs(e.normalized_error) <= k ? "1" : "0") : ""});
        groups[{r.at("split"), r.at("origin"), r.at("group"), r.at("quantity"), r.at("weighting"),
                fmt(num(r, "frequency_hz")), fmt(num(r, "bandwidth_hz")), r.at("signal"), r.at("surface"),
                fmt(std::floor(num(r, "wind_mps"))), fmt(30 * std::floor(num(r, "direction_deg") / 30))}]
            .push_back(e);
        if (r.at("split") == "validation" && r.at("origin") != "synthetic")
            ++external;
    }
    for (const auto &[key, values] : groups) {
        const auto s = error_statistics(values, k);
        auto row = key;
        for (const auto &v :
             {std::to_string(s.count), fmt(s.bias_db), fmt(s.mae_db), fmt(s.rmse_db), fmt(s.max_abs_db),
              std::to_string(s.normalized_count), s.normalized_count ? fmt(s.normalized_rmse) : "",
              std::to_string(s.within_k_count)})
            row.push_back(v);
        summaries.push_back(std::move(row));
    }
    fresh(output);
    std::filesystem::copy_file(observed, output / "observations.csv");
    std::filesystem::copy_file(predicted, output / "predictions.csv");
    write(output / "residuals.csv",
          {"id", "group", "split", "origin", "predicted_minus_measured_db", "combined_standard_u_db",
           "normalized_error", "within_k"},
          pairs);
    write(output / "groups.csv",
          {"split", "origin", "group", "quantity", "weighting", "frequency_hz", "bandwidth_hz", "signal",
           "surface", "wind_bin_lower_mps", "direction_sector_lower_deg", "count", "bias_db", "mae_db",
           "rmse_db", "max_abs_db", "normalized_count", "normalized_rmse", "within_k_count"},
          summaries);
    finish(output, "paired_errors_independent_model_measurement_uncertainties", k, external);
}
void budget_files(const std::filesystem::path &perturbations, const std::filesystem::path &correlations,
                  const std::filesystem::path &output, double k) {
    const auto input = read(perturbations, {"output_id", "parameter", "unit", "x_minus", "x0", "x_plus",
                                            "u_x", "y_minus_db", "y0_db", "y_plus_db", "provenance"});
    const auto corr = read(correlations, {"output_id", "parameter_a", "parameter_b", "rho"}, true);
    std::map<std::string, std::vector<Perturbation>> groups;
    for (const auto &r : input) {
        require(!r.at("output_id").empty(), "Missing output ID");
        groups[r.at("output_id")].push_back({r.at("parameter"), r.at("unit"), r.at("provenance"),
                                             num(r, "x_minus"), num(r, "x0"), num(r, "x_plus"), num(r, "u_x"),
                                             num(r, "y_minus_db"), num(r, "y0_db"), num(r, "y_plus_db")});
    }
    for (const auto &r : corr)
        require(groups.count(r.at("output_id")), "Unknown correlation output ID");
    std::vector<std::vector<std::string>> budgets, sensitivities;
    for (const auto &[id, p] : groups) {
        require(p.size() <= 256, "Too many uncertain parameters");
        std::map<std::string, std::size_t> index;
        for (std::size_t i = 0; i < p.size(); ++i)
            require(index.emplace(p[i].parameter, i).second, "Duplicate perturbation parameter");
        std::vector<std::vector<double>> matrix(p.size(), std::vector<double>(p.size()));
        for (std::size_t i = 0; i < p.size(); ++i)
            matrix[i][i] = 1;
        std::set<std::pair<std::size_t, std::size_t>> seen;
        for (const auto &r : corr)
            if (r.at("output_id") == id) {
                require(index.count(r.at("parameter_a")) && index.count(r.at("parameter_b")),
                        "Unknown correlation parameter");
                const auto a = index.at(r.at("parameter_a")), b = index.at(r.at("parameter_b"));
                require(a != b && seen.emplace(std::min(a, b), std::max(a, b)).second,
                        "Duplicate/diagonal correlation row");
                matrix[a][b] = matrix[b][a] = num(r, "rho");
            }
        const auto b = uncertainty_budget(p, matrix, k);
        budgets.push_back({id, fmt(b.nominal_db), fmt(b.standard_u_db), fmt(k), fmt(b.expanded_u_db),
                           fmt(b.nominal_db - b.expanded_u_db), fmt(b.nominal_db + b.expanded_u_db)});
        for (std::size_t i = 0; i < p.size(); ++i) {
            const auto &s = b.sensitivity[i];
            sensitivities.push_back({id, p[i].parameter, p[i].unit, fmt(s.derivative), fmt(s.left_derivative),
                                     fmt(s.right_derivative), fmt(s.curvature), fmt(s.contribution_u_db)});
        }
    }
    fresh(output);
    std::filesystem::copy_file(perturbations, output / "perturbations.csv");
    std::filesystem::copy_file(correlations, output / "correlations.csv");
    write(output / "budgets.csv",
          {"output_id", "nominal_db", "standard_u_db", "k", "expanded_u_db", "lower_db", "upper_db"},
          budgets);
    write(output / "sensitivities.csv",
          {"output_id", "parameter", "unit", "derivative_db_per_unit", "left_derivative", "right_derivative",
           "curvature_db_per_unit2", "individual_u_db"},
          sensitivities);
    finish(output, "first_order_correlated_uncertainty_symmetric_finite_difference", k, 0);
}
void ensemble_file(const std::filesystem::path &samples, const std::filesystem::path &output) {
    const auto input = read(samples, {"output_id", "sample_id", "level_db", "provenance"});
    std::map<std::string, std::map<std::string, double>> groups;
    for (const auto &r : input) {
        require(!r.at("output_id").empty() && !r.at("sample_id").empty() && !r.at("provenance").empty(),
                "Missing ensemble metadata");
        require(groups[r.at("output_id")].emplace(r.at("sample_id"), num(r, "level_db")).second,
                "Duplicate ensemble sample");
    }
    std::set<std::string> expected;
    for (const auto &[id, value] : groups.begin()->second)
        expected.insert(id);
    std::vector<std::vector<std::string>> result;
    for (const auto &[id, draws] : groups) {
        std::set<std::string> present;
        std::vector<double> values;
        for (const auto &[name, value] : draws) {
            present.insert(name);
            values.push_back(value);
        }
        require(present == expected,
                "Incomplete ensemble outputs; failed draws must not be silently dropped");
        const auto s = ensemble_statistics(values);
        result.push_back({id, std::to_string(s.count), fmt(s.mean_db), fmt(s.sample_sd_db), fmt(s.p025_db),
                          fmt(s.p50_db), fmt(s.p975_db), fmt(s.mean_energy_level_db)});
    }
    fresh(output);
    std::filesystem::copy_file(samples, output / "samples.csv");
    write(output / "ensemble.csv",
          {"output_id", "count", "mean_db", "sample_sd_db", "p025_db", "p50_db", "p975_db",
           "mean_energy_level_db"},
          result);
    finish(output, "equal_weight_empirical_draws_linear_quantiles_not_mean_confidence_interval", 0, 0);
}
void import_receiver_map(const std::filesystem::path &map, const std::filesystem::path &context,
                         const std::filesystem::path &output, const std::string &provenance) {
    require(!std::filesystem::exists(output), "Import output must be new");
    require(!provenance.empty(), "Missing simulation provenance");
    const auto inputs = read(map,
                             {"observer", "x_m", "y_m", "z_m", "start_receiver_s", "end_receiver_s",
                              "LAeq_turbine_dB", "LAeq_with_background_dB"},
                             false, true);
    auto rows = read(context, schema);
    std::map<std::string, const Row *> observers;
    for (const auto &r : inputs)
        require(observers.emplace(r.at("observer"), &r).second, "Duplicate receiver in map");
    std::set<std::string> ids;
    std::vector<std::vector<std::string>> result;
    for (auto &r : rows) {
        validate_record(r, true);
        require(ids.insert(r.at("id")).second, "Duplicate context ID");
        require(r.at("quantity") == "LAeq" && r.at("weighting") == "A",
                "Map import supports A-weighted LAeq only");
        require(observers.count(r.at("receiver")), "Receiver missing from map");
        const auto &m = *observers.at(r.at("receiver"));
        for (const auto &axis : {"x_m", "y_m", "z_m"})
            require(near(num(r, axis), num(m, axis)), "Receiver coordinate mismatch");
        require(near(num(r, "start_s"), num(m, "start_receiver_s")) &&
                    near(num(r, "end_s"), num(m, "end_receiver_s")),
                "Receiver averaging window mismatch");
        r["level_db"] =
            fmt(num(m, r.at("signal") == "turbine_only" ? "LAeq_turbine_dB" : "LAeq_with_background_dB"));
        r["provenance"] = provenance;
        std::vector<std::string> cells;
        for (const auto &key : schema)
            cells.push_back(r.at(key));
        result.push_back(std::move(cells));
    }
    write(output, schema, result);
}
} // namespace turbine::validation
