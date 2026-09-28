#include "trailing_edge.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
namespace aeroacoustics {
const char *tno_edge_velocity_name(TnoEdgeVelocity mode) {
    switch (mode) {
    case TnoEdgeVelocity::reference:
        return "reference";
    case TnoEdgeVelocity::input:
        return "input";
    }
    throw std::invalid_argument("Unknown TNO edge-velocity mode");
}
void TrailingEdgeOptions::validate() const {
    if (provenance.empty())
        throw std::invalid_argument("Trailing-edge configuration needs provenance");
    Parameters p;
    p.tbltemod = model;
    p.tno_edge_velocity = edge_velocity;
    p.howe = howe;
    aeroacoustics::validate(p);
    if (model != 2 && edge_velocity != TnoEdgeVelocity::reference)
        throw std::invalid_argument("TNO edge mode applies only to the TNO model");
}
void TrailingEdgeOptions::apply(Parameters &p) const {
    validate();
    p.tbltemod = model;
    p.tno_edge_velocity = edge_velocity;
    p.howe = howe;
    if (model >= 2)
        p.x_blmethod = 2;
    aeroacoustics::validate(p);
}
TrailingEdgeOptions TrailingEdgeOptions::read(const std::string &path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::invalid_argument("Cannot open trailing-edge configuration: " + path);
    std::map<std::string, std::string> fields;
    std::string line;
    while (std::getline(stream, line)) {
        if (fields.empty() && line.compare(0, 3, "\xef\xbb\xbf") == 0)
            line.erase(0, 3);
        std::istringstream row(line);
        row >> std::ws;
        if (row.eof() || row.peek() == '!')
            continue;
        std::string value, key;
        if (!(row >> std::quoted(value) >> key))
            throw std::invalid_argument("Invalid trailing-edge configuration row");
        for (auto &c : key)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        row >> std::ws;
        if (!row.eof() && row.peek() != '!')
            throw std::invalid_argument("Unexpected trailing-edge row suffix");
        if (!fields.emplace(key, value).second)
            throw std::invalid_argument("Duplicate trailing-edge field: " + key);
    }
    if (!stream.eof())
        throw std::runtime_error("Trailing-edge configuration read failed");
    auto required = [&](const char *key) -> const std::string & {
        const auto i = fields.find(key);
        if (i == fields.end())
            throw std::invalid_argument(std::string("Missing trailing-edge field: ") + key);
        return i->second;
    };
    auto number = [&](const char *key) {
        const auto &v = required(key);
        std::size_t used = 0;
        const double x = std::stod(v, &used);
        if (used != v.size() || !std::isfinite(x))
            throw std::invalid_argument("Invalid trailing-edge number");
        return x;
    };
    TrailingEdgeOptions o;
    o.provenance = required("provenance");
    const auto &m = required("model");
    if (m == "off")
        o.model = 0;
    else if (m == "bpm")
        o.model = 1;
    else if (m == "tno")
        o.model = 2;
    else if (m == "howe-chase")
        o.model = 3;
    else
        throw std::invalid_argument("Unknown trailing-edge model");
    if (o.model == 2) {
        const auto &e = required("tnoedgevelocity");
        if (e == "reference")
            o.edge_velocity = TnoEdgeVelocity::reference;
        else if (e == "input")
            o.edge_velocity = TnoEdgeVelocity::input;
        else
            throw std::invalid_argument("TNOEdgeVelocity must be reference or input");
    }
    if (o.model == 3) {
        o.howe.half_height = number("halfheightm");
        o.howe.wavelength = number("wavelengthm");
        o.howe.convection_ratio = number("convectionratio");
        o.howe.friction_ratio = number("frictionratio");
    }
    for (const auto &[key, value] : fields) {
        const bool common = key == "model" || key == "provenance";
        const bool tno = o.model == 2 && key == "tnoedgevelocity";
        const bool howe = o.model == 3 && (key == "halfheightm" || key == "wavelengthm" ||
                                           key == "convectionratio" || key == "frictionratio");
        if (!common && !tno && !howe)
            throw std::invalid_argument("Unknown/inapplicable trailing-edge field: " + key);
    }
    o.validate();
    return o;
}
} // namespace aeroacoustics
