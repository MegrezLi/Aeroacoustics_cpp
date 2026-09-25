#pragma once
#include <filesystem>
#include <string>
#include <vector>

namespace turbine::validation {
struct Residual {
    double error_db, combined_u_db, normalized_error;
    bool normalized_available;
};
struct ErrorStatistics {
    std::size_t count = 0, normalized_count = 0, within_k_count = 0;
    double bias_db = 0, mae_db = 0, rmse_db = 0, max_abs_db = 0, normalized_rmse = 0;
};
// Uncertainties here are standard uncertainties in dB. Measurement and model errors
// are assumed independent; these scores do not establish statistical independence of rows.
Residual residual(double predicted_db, double measured_db, double model_u_db, double measurement_u_db);
ErrorStatistics error_statistics(const std::vector<Residual> &, double coverage_factor);
struct Perturbation {
    std::string parameter, unit, provenance;
    double minus, nominal, plus, standard_u;
    double y_minus_db, y_nominal_db, y_plus_db;
};
struct Sensitivity {
    double derivative, left_derivative, right_derivative, curvature, contribution_u_db;
};
struct Budget {
    std::vector<Sensitivity> sensitivity;
    double nominal_db, standard_u_db, expanded_u_db;
};
// Symmetric, one-at-a-time perturbations, with a PSD correlation matrix in parameter order.
// The first-order budget includes signed covariance terms; k is not a confidence probability.
Budget uncertainty_budget(const std::vector<Perturbation> &, const std::vector<std::vector<double>> &,
                          double coverage_factor);
struct EnsembleStatistics {
    std::size_t count;
    double mean_db, sample_sd_db, p025_db, p50_db, p975_db, mean_energy_level_db;
};
// Equal-weight draws supplied by the caller, not a confidence interval for the mean.
EnsembleStatistics ensemble_statistics(std::vector<double>);

void compare_files(const std::filesystem::path &observed, const std::filesystem::path &predicted,
                   const std::filesystem::path &output, double coverage_factor);
void budget_files(const std::filesystem::path &perturbations, const std::filesystem::path &correlations,
                  const std::filesystem::path &output, double coverage_factor);
void ensemble_file(const std::filesystem::path &samples, const std::filesystem::path &output);
void import_receiver_map(const std::filesystem::path &map, const std::filesystem::path &context,
                         const std::filesystem::path &output_csv, const std::string &provenance);
} // namespace turbine::validation
