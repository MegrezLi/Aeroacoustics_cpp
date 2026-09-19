#pragma once
#include "propagation.hpp"
#include "turbine/input.hpp"
namespace turbine {
aeroacoustics::PropagationOptions read_propagation(const std::filesystem::path &);
}
