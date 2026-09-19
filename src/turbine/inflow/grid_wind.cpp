#include "turbine/wind.hpp"
#include <fstream>
#include <limits>
namespace turbine {
GridWind::GridWind(std::array<std::vector<double>, 4> axes, std::vector<Vec3> values)
    : axes_(std::move(axes)), values_(std::move(values)) {
    std::size_t count = 1;
    for (const auto &axis : axes_) {
        if (axis.empty() || axis.size() > 100000000 / count)
            throw std::invalid_argument("Invalid/oversized wind grid");
        count *= axis.size();
        for (std::size_t i = 0; i < axis.size(); ++i)
            if (!std::isfinite(axis[i]) || (i && axis[i] <= axis[i - 1]))
                throw std::invalid_argument("Wind axes must be finite and increasing");
    }
    if (values_.size() != count)
        throw std::invalid_argument("Wind grid shape mismatch");
    for (const auto &v : values_)
        for (double x : v)
            if (!std::isfinite(x))
                throw std::invalid_argument("Non-finite wind velocity");
}
std::shared_ptr<const GridWind> GridWind::read(const std::filesystem::path &path) {
    std::ifstream f(path);
    std::string magic;
    int version;
    std::array<std::size_t, 4> sizes{};
    if (!(f >> magic >> version) || magic != "AA_WIND_GRID" || version != 1)
        throw std::runtime_error("Expected AA_WIND_GRID 1: " + path.string());
    std::size_t count = 1;
    for (auto &n : sizes) {
        if (!(f >> n) || n == 0 || n > 100000000 / count)
            throw std::runtime_error("Invalid wind grid dimensions");
        count *= n;
    }
    std::array<std::vector<double>, 4> axes;
    for (std::size_t d = 0; d < 4; ++d) {
        axes[d].resize(sizes[d]);
        for (auto &v : axes[d])
            if (!(f >> v))
                throw std::runtime_error("Truncated wind grid axes");
    }
    std::vector<Vec3> values(count);
    for (auto &v : values)
        for (auto &x : v)
            if (!(f >> x))
                throw std::runtime_error("Truncated wind grid values");
    std::string extra;
    if (f >> extra)
        throw std::runtime_error("Unexpected data after wind grid");
    return std::make_shared<const GridWind>(std::move(axes), std::move(values));
}
Vec3 GridWind::at(double t, const Vec3 &p) const {
    const std::array<double, 4> query{t, p[0], p[1], p[2]};
    std::array<std::size_t, 4> lo{};
    std::array<double, 4> w{};
    for (std::size_t d = 0; d < 4; ++d) {
        const auto &a = axes_[d];
        if (!std::isfinite(query[d]))
            throw std::invalid_argument("Non-finite wind query");
        // A singleton explicitly represents a uniform dimension, including steady time.
        if (a.size() == 1)
            continue;
        if (query[d] < a.front() || query[d] > a.back())
            throw std::out_of_range("Wind query outside grid axis " + std::to_string(d));
        auto upper = std::upper_bound(a.begin(), a.end(), query[d]);
        lo[d] = upper == a.end() ? a.size() - 2 : std::size_t(upper - a.begin()) - 1;
        w[d] = (query[d] - a[lo[d]]) / (a[lo[d] + 1] - a[lo[d]]);
    }
    Vec3 result{};
    for (unsigned corner = 0; corner < 16; ++corner) {
        double weight = 1;
        std::size_t index = 0;
        for (unsigned d = 0; d < 4; ++d) {
            const bool high = (corner >> d) & 1;
            if (axes_[d].size() == 1 && high) {
                weight = 0;
                break;
            }
            weight *= high ? w[d] : 1 - w[d];
            index = index * axes_[d].size() + lo[d] + (high ? 1 : 0);
        }
        if (weight)
            result = result + weight * values_[index];
    }
    return result;
}
} // namespace turbine
