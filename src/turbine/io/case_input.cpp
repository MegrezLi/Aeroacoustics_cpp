#include "turbine/input.hpp"
#include "lookup_diagnostics.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace turbine {
namespace {
std::string lower(std::string s) {
    for (auto &c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
bool numeric(std::string s, double &x) {
    std::replace(s.begin(), s.end(), 'D', 'E');
    std::replace(s.begin(), s.end(), 'd', 'e');
    try {
        std::size_t n = 0;
        x = std::stod(s, &n);
        return n == s.size() && std::isfinite(x);
    } catch (...) {
        return false;
    }
}
std::vector<std::string> tokenize(const std::string &s) {
    std::vector<std::string> result;
    std::string word;
    char quote = 0;
    for (char c : s) {
        if (quote) {
            if (c == quote)
                quote = 0;
            else
                word += c;
        } else if ((c == '"' || c == '\'') && (word.empty() || word == "@"))
            quote = c;
        else if (c == '!')
            break;
        else if (std::isspace(static_cast<unsigned char>(c)) || c == ',') {
            if (word == "-" && result.size() >= 2)
                break;
            if (!word.empty()) {
                result.push_back(word);
                word.clear();
            }
        } else
            word += c;
    }
    if (quote)
        throw std::runtime_error("Unterminated quote in input: " + s);
    if (!word.empty())
        result.push_back(word);
    return result;
}
void require(bool value, const std::string &message) {
    if (!value)
        throw std::runtime_error(message);
}
} // namespace
InputFile::InputFile(std::filesystem::path filename) : path(std::filesystem::absolute(std::move(filename))) {
    std::ifstream stream(path);
    require(bool(stream), "Cannot open input: " + path.string());
    std::string s;
    while (std::getline(stream, s)) {
        if (rows.empty() && s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xef)
            s.erase(0, 3);
        auto tokens = tokenize(s);
        if (tokens.size() >= 2 && tokens[0].find_first_not_of("-=#") != std::string::npos &&
            std::isalpha(static_cast<unsigned char>(tokens[1][0]))) {
            const auto key = lower(tokens[1]);
            if (labels_.count(key) == 0)
                labels_.emplace(key, rows.size());
        }
        rows.push_back(std::move(tokens));
    }
}
bool InputFile::has(const std::string &key) const { return labels_.count(lower(key)) != 0; }
std::size_t InputFile::line(const std::string &key) const {
    const auto it = labels_.find(lower(key));
    if (it == labels_.end())
        throw std::runtime_error(path.string() + ": missing field " + key);
    return it->second;
}
std::string InputFile::value(const std::string &key) const { return rows.at(line(key)).at(0); }
double InputFile::number(const std::string &key) const {
    double x = 0;
    require(numeric(value(key), x), path.string() + ": invalid number for " + key + ": " + value(key));
    return x;
}
double InputFile::number(const std::string &key, double fallback) const {
    return !has(key) || lower(value(key)) == "default" ? fallback : number(key);
}
int InputFile::integer(const std::string &key) const {
    double x = number(key);
    require(x == std::floor(x) && std::abs(x) < 1e8, "Invalid integer: " + key);
    return static_cast<int>(x);
}
bool InputFile::flag(const std::string &key) const {
    auto s = lower(value(key));
    if (s == "true" || s == "t")
        return true;
    if (s == "false" || s == "f")
        return false;
    throw std::runtime_error(path.string() + ": invalid logical field " + key);
}
std::filesystem::path InputFile::file(const std::string &key) const {
    return (path.parent_path() / value(key)).lexically_normal();
}
std::vector<std::vector<double>> InputFile::table_after(const std::string &key, std::size_t count,
                                                        std::size_t columns) const {
    std::vector<std::vector<double>> result;
    for (std::size_t i = line(key) + 1; i < rows.size() && result.size() < count; ++i) {
        const auto &row = rows[i];
        double first = 0;
        if (row.empty() || !numeric(row[0], first))
            continue;
        require(row.size() >= columns, path.string() + ": truncated table after " + key);
        std::vector<double> values;
        for (std::size_t c = 0; c < columns; ++c) {
            double x;
            require(numeric(row[c], x), path.string() + ": invalid table value after " + key);
            values.push_back(x);
        }
        result.push_back(std::move(values));
    }
    require(result.size() == count, path.string() + ": incomplete table after " + key);
    return result;
}
std::vector<std::filesystem::path> InputFile::files_after(const std::string &key, std::size_t count) const {
    std::vector<std::filesystem::path> result;
    for (std::size_t i = line(key); i < rows.size() && result.size() < count; ++i) {
        if (rows[i].empty())
            continue;
        result.push_back((path.parent_path() / rows[i][0]).lexically_normal());
    }
    require(result.size() == count, "Incomplete filename list: " + key);
    return result;
}
Airfoil::Airfoil(const std::filesystem::path &filename) : input(filename) {
    require(input.integer("NumTabs") == 1, "Only one airfoil table is supported");
    // AirfoilInfo.f90 DefaultInterpOrd=1 in the pinned version, not 3 as the file comment claims.
    require(input.number("InterpOrd", 1) == 1, "Cubic airfoil interpolation is not implemented yet");
    for (const auto &row : input.table_after("NumAlf", input.integer("NumAlf"), 4)) {
        require(alpha.empty() || row[0] * deg > alpha.back(), "Airfoil angles must increase");
        alpha.push_back(row[0] * deg);
        coefficients.push_back({row[1], row[2], row[3]});
    }
    require(alpha.size() >= 2, "An airfoil table needs at least two angles");
    auto coord = input.value("NumCoords");
    if (!coord.empty() && coord[0] == '@') {
        std::ifstream f(input.path.parent_path() / coord.substr(1));
        require(bool(f), "Cannot open airfoil coordinates: " + coord);
        std::string row;
        std::vector<std::vector<double>> numbers;
        while (std::getline(f, row)) {
            auto words = tokenize(row);
            if (words.empty())
                continue;
            std::vector<double> a;
            double x;
            for (const auto &w : words) {
                if (!numeric(w, x))
                    break;
                a.push_back(x);
            }
            if (!a.empty())
                numbers.push_back(a);
        }
        // External coordinate format: count; reference point; coordinate pairs.
        require(numbers.size() >= 3 && numbers[0].size() == 1 && numbers[1].size() == 2,
                "Invalid airfoil coordinates: " + coord);
        const auto n = static_cast<std::size_t>(numbers[0][0]);
        reference = {numbers[1][0], numbers[1][1]};
        require(numbers.size() == n + 1, "Airfoil coordinate count mismatch: " + coord);
        for (std::size_t i = 2; i < numbers.size(); ++i) {
            require(numbers[i].size() == 2, "Invalid coordinate pair");
            coordinates.push_back({numbers[i][0], numbers[i][1]});
        }
    } else
        require(coord == "0", "Inline airfoil coordinates are not implemented yet");
}
Coefficients Airfoil::at(double a) const {
    if (!std::isfinite(a))
        throw std::invalid_argument("Non-finite airfoil angle: " + input.path.string());
    a = std::remainder(a, 2 * pi);
    // Convert only exceptional queries to degrees; valid hot-path queries do not allocate.
    if (a < alpha.front() || a > alpha.back())
        diagnostics::check_lookup(input.path.string(), "alpha_deg", a / deg, alpha.front() / deg, alpha.back() / deg);
    if (a <= alpha.front())
        return coefficients.front();
    if (a >= alpha.back())
        return coefficients.back();
    const auto j = static_cast<std::size_t>(std::upper_bound(alpha.begin(), alpha.end(), a) - alpha.begin());
    const double w = (a - alpha[j - 1]) / (alpha[j] - alpha[j - 1]);
    const auto &l = coefficients[j - 1];
    const auto &r = coefficients[j];
    return {l.cl + w * (r.cl - l.cl), l.cd + w * (r.cd - l.cd), l.cm + w * (r.cm - l.cm)};
}
SteadyWind::SteadyWind(const InputFile &f)
    : speed(f.number("HWindSpeed")), reference_height(f.number("RefHt")), exponent(f.number("PLExp")),
      propagation(f.number("PropagationDir") * deg), upflow(f.number("VFlowAng") * deg) {
    require(f.integer("WindType") == 1, "Only steady InflowWind WindType=1 is supported");
    require(reference_height > 0 && speed >= 0, "Invalid steady wind input");
}
Vec3 SteadyWind::at(const Vec3 &p) const {
    // Avoid constructing require()'s std::string on every valid node query.
    if (!(p[2] > 0 || exponent == 0))
        throw std::runtime_error("Power law wind requested below ground");
    double v = exponent == 0 ? speed : speed * std::pow(p[2] / reference_height, exponent);
    return {v * std::cos(upflow) * std::cos(propagation), -v * std::cos(upflow) * std::sin(propagation),
            v * std::sin(upflow)};
}
Case::Case(const std::filesystem::path &filename)
    : primary(filename), structure(primary.file("EDFile")), aero(primary.file("AeroFile")),
      inflow(primary.file("InflowFile")), acoustic(aero.file("AA_InputFile")),
      blade_structure(structure.file("BldFile(1)")), wind(inflow), dt(primary.number("DT")),
      duration(primary.number("TMax")), rho(primary.number("AirDens")), nu(primary.number("KinVisc")),
      sound_speed(primary.number("SpdSound")), gravity(primary.number("Gravity")) {
    validate_scope();
    rho = aero.number("AirDens", rho);
    nu = aero.number("KinVisc", nu);
    sound_speed = aero.number("SpdSound", sound_speed);
    require(rho > 0 && nu > 0 && sound_speed > 0, "Invalid aerodynamic environment");
    InputFile blade(aero.file("ADBlFile(1)"));
    for (const auto &row : blade.table_after("NumBlNds", blade.integer("NumBlNds"), 7)) {
        require(stations.empty() || row[0] > stations.back().span, "Blade span must increase");
        require(row[5] > 0 && row[6] == std::floor(row[6]), "Invalid blade chord or airfoil ID");
        stations.push_back(
            {row[0], row[1], row[2], row[3] * deg, row[4] * deg, row[5], static_cast<int>(row[6]) - 1});
    }
    for (const auto &f : aero.files_after("AFNames", aero.integer("NumAFfiles")))
        airfoils.emplace_back(f);
    for (const auto &s : stations)
        require(s.airfoil >= 0 && static_cast<std::size_t>(s.airfoil) < airfoils.size(),
                "Airfoil index outside table list");
}
void Case::validate_scope() const {
    require(primary.integer("CompElast") == 1 && primary.integer("CompInflow") == 1 &&
                primary.integer("CompAero") == 2,
            "Unsupported module configuration");
    for (auto key :
         {"CompServo", "CompSeaSt", "CompHydro", "CompSub", "CompMooring", "CompIce", "CompSoil", "MHK"})
        require(primary.integer(key) == 0, std::string("Unsupported enabled module: ") + key);
    require(primary.integer("NRotors") == 1 && !primary.flag("MirrorRotor"),
            "Only one normal rotor is supported");
    require(structure.integer("NumBl") == 3, "Only three blades are supported");
    for (auto key : {"TTDspFA", "TTDspSS", "PtfmSurge", "PtfmSway", "PtfmHeave", "PtfmRoll", "PtfmPitch",
                     "PtfmYaw", "PtfmRefxt", "PtfmRefyt", "PtfmRefzt"})
        require(structure.number(key) == 0, std::string("Unsupported fixed structural offset: ") + key);
    for (int b = 1; b <= 3; ++b) {
        const auto suffix = "(" + std::to_string(b) + ")";
        require(structure.number("TipMass" + suffix) == 0, "Nonzero tip-brake masses are not supported");
        require(structure.file("BldFile" + suffix) == structure.file("BldFile(1)"),
                "All blades must use the same structural property file");
        require(aero.file("ADBlFile" + suffix) == aero.file("ADBlFile(1)"),
                "All blades must use the same aerodynamic station file");
    }
    for (auto key : {"FlapDOF1", "FlapDOF2", "EdgeDOF"})
        require(structure.flag(key), std::string("Expected active blade DOF: ") + key);
    for (auto key :
         {"PitchDOF", "TeetDOF", "DrTrDOF", "GenDOF", "YawDOF", "TwFADOF1", "TwFADOF2", "TwSSDOF1",
          "TwSSDOF2", "PtfmSgDOF", "PtfmSwDOF", "PtfmHvDOF", "PtfmRDOF", "PtfmPDOF", "PtfmYDOF", "Furling"})
        require(!structure.flag(key), std::string("Unsupported active DOF: ") + key);
    require(aero.integer("Wake_Mod") == 1 && aero.integer("BEM_Mod") == 1 && aero.integer("DBEMT_Mod") == 0,
            "Unsupported wake model");
    require(aero.integer("UA_Mod") == 3 && aero.flag("FLookup") && aero.integer("AFTabMod") == 1,
            "Unsupported unsteady airfoil model");
    require(aero.integer("TwrPotent") == 0 && aero.integer("TwrShadow") == 0 && !aero.flag("TwrAero") &&
                !aero.flag("NacelleDrag") && !aero.flag("TFinAero"),
            "Unsupported tower, nacelle or tail aerodynamics");
    require(!aero.flag("SectAvg") && !aero.flag("SkewMomCorr"),
            "Unsupported BEM averaging or skew momentum correction");
    require(aero.number("DTAero", dt) == dt, "Aerodynamic sub-stepping is not supported");
    require(aero.integer("Skew_Mod") >= 0 && aero.integer("Skew_Mod") <= 1 &&
                (aero.number("SkewRedistr_Mod", 1) == 0 || aero.number("SkewRedistr_Mod", 1) == 1),
            "Unsupported skew model");
    require(aero.number("UAStartRad", 0) == 0 && aero.number("UAEndRad", 1) == 1,
            "UA radial clipping is not supported");
    require(aero.integer("InCol_Alfa") == 1 && aero.integer("InCol_Cl") == 2 &&
                aero.integer("InCol_Cd") == 3 && aero.integer("InCol_Cm") == 4,
            "Unsupported airfoil table columns");
    require(primary.number("RhoInf") >= 0 && primary.number("RhoInf") <= 1, "RhoInf must be in [0,1]");
    require(!primary.flag("Linearize") && !primary.flag("CalcSteady"),
            "Only time-domain simulation is supported");
    require(aero.flag("CompAA"), "CompAA must be enabled for the acoustic turbine solver");
    require(acoustic.integer("NrOutFile") >= 1 && acoustic.integer("NrOutFile") <= 4,
            "NrOutFile must be in [1,4]");
    require(acoustic.integer("TICalcMeth") == 1,
            "The coupled turbine solver currently supports TICalcMeth=1");
    const double acoustic_dt = acoustic.number("DT_AA", dt);
    require(acoustic_dt >= dt && std::abs(acoustic_dt / dt - std::round(acoustic_dt / dt)) < 1e-9,
            "DT_AA must be a positive multiple of DT");
    require(dt > 0 && duration >= 0 && rho > 0 && nu > 0 && sound_speed > 0,
            "Invalid time or environment inputs");
}
} // namespace turbine
