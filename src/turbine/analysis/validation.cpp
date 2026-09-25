#include "turbine/validation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace turbine::validation {
namespace {
void require(bool ok, const char *message) {
    if (!ok)
        throw std::invalid_argument(message);
}
void finite(double v) { require(std::isfinite(v), "Non-finite validation input/result"); }
void factor(double k) {
    finite(k);
    require(k > 0, "Coverage factor must be positive");
}
double checked(long double v) {
    double d = static_cast<double>(v);
    finite(d);
    return d;
}
} // namespace
Residual residual(double predicted, double measured, double model_u, double measurement_u) {
    for (double v : {predicted, measured, model_u, measurement_u})
        finite(v);
    require(model_u >= 0 && measurement_u >= 0, "Standard uncertainty must be nonnegative");
    const double error = checked(static_cast<long double>(predicted) - measured);
    const double u = std::hypot(model_u, measurement_u);
    finite(u);
    const double normalized = u > 0 ? error / u : 0;
    finite(normalized);
    return {error, u, normalized, u > 0};
}
ErrorStatistics error_statistics(const std::vector<Residual> &values, double k) {
    factor(k);
    require(!values.empty(), "Empty validation group");
    ErrorStatistics out;
    long double sum = 0, absolute = 0, square = 0, normalized_square = 0;
    for (const auto &v : values) {
        finite(v.error_db);
        finite(v.combined_u_db);
        require(v.combined_u_db >= 0 && v.normalized_available == (v.combined_u_db > 0),
                "Invalid residual uncertainty");
        const long double e = v.error_db;
        sum += e;
        absolute += std::abs(e);
        square += e * e;
        out.max_abs_db = std::max(out.max_abs_db, std::abs(v.error_db));
        if (v.normalized_available) {
            // Recompute rather than trust a caller-supplied cached normalized value.
            const long double z = e / v.combined_u_db;
            normalized_square += z * z;
            ++out.normalized_count;
            if (std::abs(z) <= k)
                ++out.within_k_count;
        }
    }
    out.count = values.size();
    out.bias_db = checked(sum / out.count);
    out.mae_db = checked(absolute / out.count);
    out.rmse_db = checked(std::sqrt(square / out.count));
    if (out.normalized_count)
        out.normalized_rmse = checked(std::sqrt(normalized_square / out.normalized_count));
    return out;
}
Budget uncertainty_budget(const std::vector<Perturbation> &p, const std::vector<std::vector<double>> &r,
                          double k) {
    factor(k);
    const auto n = p.size();
    require(n > 0 && n <= 256 && r.size() == n, "Invalid uncertainty matrix size (1..256)");
    for (const auto &row : r)
        require(row.size() == n, "Non-square correlation matrix");
    std::vector<std::vector<long double>> l(n, std::vector<long double>(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            finite(r[i][j]);
            require(std::abs(r[i][j]) <= 1 && std::abs(r[i][j] - r[j][i]) <= 1e-12,
                    "Invalid correlation symmetry/range");
        }
        require(std::abs(r[i][i] - 1) <= 1e-12, "Correlation diagonal must equal one");
        for (std::size_t j = 0; j <= i; ++j) {
            long double v = r[i][j];
            for (std::size_t h = 0; h < j; ++h)
                v -= l[i][h] * l[j][h];
            if (i == j) {
                require(v >= -1e-12L, "Correlation matrix must be positive semidefinite");
                l[i][j] = std::sqrt(std::max(0.L, v));
            } else if (l[j][j] > 1e-14L)
                l[i][j] = v / l[j][j];
            else
                require(std::abs(v) <= 1e-12L, "Inconsistent singular correlation matrix");
        }
    }
    Budget out;
    out.nominal_db = p.front().y_nominal_db;
    std::set<std::string> names;
    std::vector<long double> signed_u;
    for (const auto &v : p) {
        require(!v.parameter.empty() && !v.unit.empty() && !v.provenance.empty() &&
                    names.insert(v.parameter).second,
                "Missing/duplicate parameter metadata");
        for (double x : {v.minus, v.nominal, v.plus, v.standard_u, v.y_minus_db, v.y_nominal_db, v.y_plus_db})
            finite(x);
        const long double a = static_cast<long double>(v.nominal) - v.minus,
                          b = static_cast<long double>(v.plus) - v.nominal;
        require(a > 0 && b > 0 && std::abs(a - b) <= 1e-8L * std::max(a, b),
                "Perturbations must bracket nominal symmetrically");
        require(v.standard_u >= 0 && std::abs(v.y_nominal_db - out.nominal_db) <= 1e-8,
                "Inconsistent baseline/uncertainty");
        const long double left = (static_cast<long double>(v.y_nominal_db) - v.y_minus_db) / a;
        const long double right = (static_cast<long double>(v.y_plus_db) - v.y_nominal_db) / b;
        const long double slope = (static_cast<long double>(v.y_plus_db) - v.y_minus_db) / (a + b);
        signed_u.push_back(slope * v.standard_u);
        out.sensitivity.push_back({checked(slope), checked(left), checked(right),
                                   checked(2 * (right - left) / (a + b)),
                                   checked(std::abs(signed_u.back()))});
    }
    // ||L^T (c*u)|| avoids a negative variance from cancellation in highly correlated budgets.
    long double variance = 0;
    for (std::size_t j = 0; j < n; ++j) {
        long double value = 0;
        for (std::size_t i = j; i < n; ++i)
            value += l[i][j] * signed_u[i];
        variance += value * value;
    }
    out.standard_u_db = checked(std::sqrt(variance));
    out.expanded_u_db = checked(k * std::sqrt(variance));
    return out;
}
EnsembleStatistics ensemble_statistics(std::vector<double> v) {
    require(v.size() >= 2, "Ensemble needs at least two equal-weight draws");
    for (double x : v)
        finite(x);
    std::sort(v.begin(), v.end());
    long double mean = 0, ss = 0, power = 0;
    std::size_t n = 0;
    for (double x : v) {
        ++n;
        long double d = static_cast<long double>(x) - mean;
        mean += d / n;
        ss += d * (x - mean);
    }
    for (double x : v)
        power += std::pow(10.L, (static_cast<long double>(x) - v.back()) / 10);
    auto quantile = [&](double p) {
        const double index = (v.size() - 1) * p;
        const auto i = std::size_t(index);
        return checked(static_cast<long double>(v[i]) +
                       (index - i) * (static_cast<long double>(v[std::min(i + 1, v.size() - 1)]) - v[i]));
    };
    return {n,
            checked(mean),
            checked(std::sqrt(ss / (n - 1))),
            quantile(.025),
            quantile(.5),
            quantile(.975),
            checked(v.back() + 10 * std::log10(power / n))};
}
} // namespace turbine::validation
