#include "turbine/validation.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace turbine::validation;
void check(bool ok) {
    if (!ok)
        throw std::runtime_error("Validation analytic check failed");
}
void close(double a, double b) { check(std::abs(a - b) < 1e-10 * std::max(1., std::abs(b))); }
template <class F> void rejects(F f) {
    bool failed = false;
    try {
        f();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed);
}
int main() {
    try {
        const auto r = residual(43, 40, 3, 4);
        close(r.error_db, 3);
        close(r.combined_u_db, 5);
        close(r.normalized_error, .6);
        auto s =
            error_statistics({residual(42, 40, 0, 1), residual(38, 40, 0, 1), residual(40, 40, 0, 0)}, 2);
        close(s.bias_db, 0);
        close(s.mae_db, 4. / 3);
        close(s.rmse_db, std::sqrt(8. / 3));
        check(s.normalized_count == 2 && s.within_k_count == 2);
        rejects([] { residual(0, 0, -1, 1); });
        rejects([] { residual(INFINITY, 0, 0, 1); });
        rejects([] { error_statistics({}, 2); });
        std::vector<Perturbation> p = {{"a", "m", "synthetic", -1, 0, 1, 2, 8, 10, 12},
                                       {"b", "s", "synthetic", -1, 0, 1, 1, 13, 10, 7}};
        auto b = uncertainty_budget(p, {{1, .5}, {.5, 1}}, 2);
        close(b.standard_u_db, std::sqrt(13));
        close(b.expanded_u_db, 2 * std::sqrt(13));
        close(b.sensitivity[0].derivative, 2);
        close(b.sensitivity[1].derivative, -3);
        close(uncertainty_budget(p, {{1, 1}, {1, 1}}, 1).standard_u_db, 1);
        close(uncertainty_budget(p, {{1, -1}, {-1, 1}}, 1).standard_u_db, 7);
        rejects([&] { uncertainty_budget(p, {{1, 2}, {2, 1}}, 2); });
        rejects([&] { uncertainty_budget(p, {{1, 0}, {1, 1}}, 2); });
        auto three = p;
        three.push_back({"c", "m", "synthetic", -1, 0, 1, 1, 9, 10, 11});
        rejects([&] { uncertainty_budget(three, {{1, -.9, -.9}, {-.9, 1, -.9}, {-.9, -.9, 1}}, 2); });
        auto bad = p;
        bad[1].y_nominal_db = 11;
        rejects([&] { uncertainty_budget(bad, {{1, 0}, {0, 1}}, 2); });
        bad = p;
        bad[0].plus = 2;
        rejects([&] { uncertainty_budget(bad, {{1, 0}, {0, 1}}, 2); });
        b = uncertainty_budget({{"quadratic", "m", "synthetic", 1, 2, 3, .2, 1, 4, 9}}, {{1}}, 2);
        close(b.sensitivity[0].derivative, 4);
        close(b.sensitivity[0].curvature, 2);
        close(b.standard_u_db, .8);
        const auto e = ensemble_statistics({0, 10, 20});
        close(e.mean_db, 10);
        close(e.sample_sd_db, 10);
        close(e.p025_db, .5);
        close(e.p50_db, 10);
        close(e.p975_db, 19.5);
        close(e.mean_energy_level_db, 10 * std::log10(37));
        close(ensemble_statistics({4000, 4000}).mean_energy_level_db, 4000);
        rejects([] { ensemble_statistics({1}); });
        rejects([] { ensemble_statistics({1, NAN}); });
        std::cout << "E8 residual, correlated/singular uncertainty, curvature and ensemble analytic checks "
                     "passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
