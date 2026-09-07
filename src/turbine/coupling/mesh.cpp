#include "turbine/mesh.hpp"
#include <limits>
namespace turbine {
MotionMap::MotionMap(std::vector<ReferenceNode> source, std::vector<ReferenceNode> destination)
    : source_(std::move(source)), destination_(std::move(destination)) {
    if (source_.size() < 2)
        throw std::runtime_error("Motion source needs two nodes");
    for (const auto &d : destination_) {
        double best = std::numeric_limits<double>::infinity();
        Projection p{};
        for (std::size_t j = 0; j + 1 < source_.size(); ++j) {
            const Vec3 edge = source_[j + 1].position - source_[j].position;
            const double w =
                std::clamp(dot(d.position - source_[j].position, edge) / dot(edge, edge), 0.0, 1.0);
            const Vec3 delta = d.position - source_[j].position - w * edge;
            const double distance = dot(delta, delta);
            if (distance < best) {
                best = distance;
                p = {j, w};
            }
        }
        map_.push_back(p);
    }
}
std::vector<Motion> MotionMap::transfer(const std::vector<Motion> &source) const {
    if (source.size() != source_.size())
        throw std::runtime_error("Motion source count mismatch");
    std::vector<Motion> result(destination_.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        auto &y = result[i];
        const auto p = map_[i];
        std::array<Matrix3, 2> orientation;
        for (int k = 0; k < 2; ++k) {
            const auto j = p.element + k;
            const double w = k ? p.weight : 1 - p.weight;
            const auto rotation = multiply(transpose(source[j].orientation), source_[j].orientation);
            y.position =
                y.position +
                w * (source[j].position + multiply(rotation, destination_[i].position - source_[j].position));
            orientation[k] = multiply(destination_[i].orientation,
                                      multiply(transpose(source_[j].orientation), source[j].orientation));
        }
        y.orientation = interpolate_rotation(orientation[0], orientation[1], p.weight);
        for (int k = 0; k < 2; ++k) {
            const auto j = p.element + k;
            const double w = k ? p.weight : 1 - p.weight;
            y.velocity = y.velocity + w * (source[j].velocity + cross(source[j].position - y.position,
                                                                      source[j].angular_velocity));
            y.angular_velocity = y.angular_velocity + w * source[j].angular_velocity;
        }
    }
    return result;
}
LoadMap::LoadMap(std::vector<Vec3> source, std::vector<Vec3> destination)
    : source_count_(source.size()), destination_count_(destination.size()) {
    if (source.size() < 2 || destination.empty())
        throw std::runtime_error("Empty load mapping mesh");
    auto nearest = [&](const Vec3 &p) {
        std::size_t best = 0;
        double distance = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < destination.size(); ++j) {
            const Vec3 d = p - destination[j];
            if (dot(d, d) < distance) {
                best = j;
                distance = dot(d, d);
            }
        }
        return best;
    };
    for (std::size_t i = 0; i + 1 < source.size(); ++i) {
        const Vec3 edge = source[i + 1] - source[i];
        const double length = norm(edge);
        std::vector<double> cuts{0, 1};
        if (length < 1e-10)
            throw std::runtime_error("Zero-length load element");
        for (const auto &d : destination) {
            const double w = dot(d - source[i], edge) / (length * length);
            if (w > 1.4901161193847656e-8 && w < 1 - 1.4901161193847656e-8)
                cuts.push_back(w);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(
            std::unique(cuts.begin(), cuts.end(), [](double a, double b) { return std::abs(a - b) < 1e-10; }),
            cuts.end());
        for (std::size_t j = 0; j + 1 < cuts.size(); ++j) {
            const double a = cuts[j], b = cuts[j + 1];
            segments_.push_back(
                {i, a, b, (b - a) * length, nearest(source[i] + a * edge), nearest(source[i] + b * edge)});
        }
    }
}
std::vector<PointLoad> LoadMap::transfer(const std::vector<PointLoad> &loads, const std::vector<Vec3> &source,
                                         const std::vector<Vec3> &destination) const {
    if (loads.size() != source_count_ || source.size() != source_count_ ||
        destination.size() != destination_count_)
        throw std::runtime_error("Load mesh count mismatch");
    std::vector<PointLoad> result(destination_count_);
    for (const auto &s : segments_) {
        const auto i = s.source;
        const Vec3 fa = (1 - s.a) * loads[i].force + s.a * loads[i + 1].force,
                   fb = (1 - s.b) * loads[i].force + s.b * loads[i + 1].force;
        const Vec3 ma = (1 - s.a) * loads[i].moment + s.a * loads[i + 1].moment,
                   mb = (1 - s.b) * loads[i].moment + s.b * loads[i + 1].moment;
        const Vec3 pa = (1 - s.a) * source[i] + s.a * source[i + 1],
                   pb = (1 - s.b) * source[i] + s.b * source[i + 1];
        const Vec3 force_a = (s.length / 6) * (2 * fa + fb), force_b = (s.length / 6) * (fa + 2 * fb),
                   couple = (s.length / 12) * cross(pb - pa, fa + fb);
        auto &a = result[s.dest_a];
        auto &b = result[s.dest_b];
        a.force = a.force + force_a;
        b.force = b.force + force_b;
        a.moment =
            a.moment + (s.length / 6) * (2 * ma + mb) + couple + cross(pa - destination[s.dest_a], force_a);
        b.moment =
            b.moment + (s.length / 6) * (ma + 2 * mb) - couple + cross(pb - destination[s.dest_b], force_b);
    }
    return result;
}
} // namespace turbine
