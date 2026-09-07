#pragma once
#include "aeroacoustics.hpp"
#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace turbine {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double deg = pi / 180.0;
using aeroacoustics::Vec3;

// OpenFAST value/label input files, including quoted paths with spaces.
// Labels are indexed only when preceded by a value; tables remain raw rows.
class InputFile {
  public:
    explicit InputFile(std::filesystem::path path);
    std::filesystem::path path;
    std::vector<std::vector<std::string>> rows;
    bool has(const std::string &key) const;
    std::size_t line(const std::string &key) const;
    std::string value(const std::string &key) const;
    double number(const std::string &key) const;
    double number(const std::string &key, double default_value) const;
    int integer(const std::string &key) const;
    bool flag(const std::string &key) const;
    std::filesystem::path file(const std::string &key) const;
    std::vector<std::vector<double>> table_after(const std::string &key, std::size_t count,
                                                 std::size_t columns) const;
    std::vector<std::filesystem::path> files_after(const std::string &key, std::size_t count) const;

  private:
    std::map<std::string, std::size_t> labels_;
};
struct Coefficients {
    double cl = 0, cd = 0, cm = 0;
};
struct Airfoil {
    InputFile input;
    std::vector<double> alpha;
    std::vector<Coefficients> coefficients;
    std::vector<std::array<double, 2>> coordinates;
    std::array<double, 2> reference{{.25, 0}};
    explicit Airfoil(const std::filesystem::path &);
    Coefficients at(double alpha_rad) const;
};
struct BladeStation {
    double span, curve, sweep, curve_angle, twist, chord;
    int airfoil;
};
struct SteadyWind {
    double speed, reference_height, exponent, propagation, upflow;
    explicit SteadyWind(const InputFile &);
    Vec3 at(const Vec3 &) const;
};
struct Case {
    InputFile primary, structure, aero, inflow, acoustic, blade_structure;
    SteadyWind wind;
    std::vector<BladeStation> stations;
    std::vector<Airfoil> airfoils;
    double dt, duration, rho, nu, sound_speed, gravity;
    explicit Case(const std::filesystem::path &);
    void validate_scope() const;
};
} // namespace turbine
