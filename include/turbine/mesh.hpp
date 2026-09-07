#pragma once
#include "turbine/structure.hpp"
namespace turbine {
struct ReferenceNode {
    Vec3 position;
    Matrix3 orientation;
};
struct Projection {
    std::size_t element;
    double weight;
};
class MotionMap {
  public:
    MotionMap(std::vector<ReferenceNode> source, std::vector<ReferenceNode> destination);
    std::vector<Motion> transfer(const std::vector<Motion> &source) const;

  private:
    std::vector<ReferenceNode> source_, destination_;
    std::vector<Projection> map_;
};
// Distributed line loads (N/m, Nm/m) to structural point loads (N, Nm).
class LoadMap {
  public:
    LoadMap(std::vector<Vec3> source, std::vector<Vec3> destination);
    std::vector<PointLoad> transfer(const std::vector<PointLoad> &distributed,
                                    const std::vector<Vec3> &source_positions,
                                    const std::vector<Vec3> &destination_positions) const;

  private:
    struct Segment {
        std::size_t source;
        double a, b, length;
        std::size_t dest_a, dest_b;
    };
    std::vector<Segment> segments_;
    std::size_t source_count_, destination_count_;
};
} // namespace turbine
