#pragma once
#include "turbine/input.hpp"
#include <memory>
namespace turbine {
// Copy/freeze once; clients cannot mutate the owned Case or invalidate its airfoils.
class TurbineModel {
    std::shared_ptr<const Case> data_;

  public:
    explicit TurbineModel(const Case &data) : data_(std::make_shared<const Case>(data)) {}
    const Case &data() const { return *data_; }
    bool same_model(const TurbineModel &other) const noexcept { return data_ == other.data_; }

  private:
    friend class Rotor;
    std::shared_ptr<const Airfoil> airfoil(std::size_t index) const {
        return std::shared_ptr<const Airfoil>(data_, &data_->airfoils.at(index));
    }
};
} // namespace turbine
