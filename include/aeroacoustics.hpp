#pragma once
#include "acoustic_quantities.hpp"
#include "mechanisms.hpp"
#include <array>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace aeroacoustics {
class OutdoorPropagation;
using Spectrum = std::vector<double>;
using Mechanisms = std::array<Spectrum, mechanism_count>;
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>; // row-major global-to-local
struct Parameters {
    Spectrum freqlist{openfast_centers_hz.begin(), openfast_centers_hz.end()};
    double spdsound = 340., kinvisc = 1.48e-5, airdens = 1.225, lturb = 40., alprat = 1., ti = .1, avgv = 8.;
    int x_blmethod = 1, itrip = 1, timod = 1, tbltemod = 1, lammod = 0, tipmod = 0, bluntmod = 0;
    bool round = true, aweighting = false;
};
enum class TrailingEdgeModel { off, bpm, tno_with_bpm_separation };
enum class InflowModel { off, lowson, lowson_guidati };
struct SourceSelection {
    TrailingEdgeModel trailing;
    InflowModel inflow;
    bool laminar, bluntness, tip;
    explicit SourceSelection(const Parameters &);
};
struct BoundaryLayer {
    std::array<double, 2> dstar{}, d99{}, cf{}, edge_velocity_ratio{{1., 1.}};
};
struct Geometry {
    double distance = 1.22, theta = 90., phi = 90.;
};
struct Section {
    double chord = .2286, speed = 63.92, alpha_deg = 3., span = .509, stall_deg = 12.5;
    Geometry trailing{};
    std::optional<Geometry> leading;
    BoundaryLayer bl{};
    double ti_section = -1., te_thickness = .001, te_angle = 14., thickness_1p = .02, thickness_10p = .12;
    bool is_tip = true;
};
class PreparedBLTable;
struct BLTable {
    Spectrum aoa, reynolds;
    std::vector<std::array<double, 8>> values; // Re-major, then AoA
    std::string source_name = "in-memory BL table";
    static BLTable read(const std::string &path);
    PreparedBLTable prepare() const;
    // Public mutable tables retain validation on every interpolation.
    BoundaryLayer interpolate(double alpha_deg, double re, double chord) const;
};
// Owns a validated copy: later changes to a BLTable cannot invalidate it.
class PreparedBLTable {
    BLTable table_;

  public:
    explicit PreparedBLTable(BLTable);
    BoundaryLayer interpolate(double alpha_deg, double re, double chord) const;
};
struct Node {
    Section section;
    Vec3 aero_center{}, inflow{};
    Mat3 global_to_local{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
    std::array<double, 2> airfoil_reference{{.25, 0}};
    std::size_t blade_number = 0, node_number = 0; // optional diagnostic IDs, one-based
};
using Snapshot = std::vector<std::vector<Mechanisms>>; // observer, node,
                                                       // mechanism, frequency
// Callback borrows one ordered observer block; do not retain it or reenter the workspace.
using ObserverBlockCallback = std::function<void(std::size_t first_observer, const Snapshot &)>;
// Owns validated options, prepared sources and reusable spectrum buffers.
class AcousticWorkspace {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    explicit AcousticWorkspace(Parameters);
    void set_propagation(std::shared_ptr<const OutdoorPropagation>);
    ~AcousticWorkspace();
    AcousticWorkspace(const AcousticWorkspace &);
    AcousticWorkspace &operator=(const AcousticWorkspace &);
    AcousticWorkspace(AcousticWorkspace &&) noexcept;
    AcousticWorkspace &operator=(AcousticWorkspace &&) noexcept;
    // View remains valid until the next evaluate() or destruction.
    const Snapshot &evaluate(const std::vector<Node> &, const std::vector<Vec3> &observers);
    void evaluate_blocks(const std::vector<Node> &, const std::vector<Vec3> &observers,
                         const ObserverBlockCallback &, std::size_t block_size = 1);
};
std::pair<std::size_t, Spectrum> blade_elements(const Spectrum &span, double percentage = 100.);
std::array<double, 2> guidati_thickness(const std::vector<std::array<double, 2>> &coords);
Snapshot snapshot_spectrum(const Parameters &, const std::vector<Node> &, const std::vector<Vec3> &);
class TurbulenceState {
    Spectrum span_, buffers_;
    std::vector<std::size_t> radial_, counts_;
    std::size_t samples_, blades_;
    double height_, ti_, avgv_;
    int method_;

  public:
    Spectrum values; // blade-major, then radial node
    TurbulenceState(Spectrum span, std::size_t blades, double dt, double hub_height, int method = 1,
                    double ti = .1, double avgv = 8.);
    void update(const Spectrum &vrel, const std::vector<Vec3> &inflow, const std::vector<Vec3> &leading);
};
class AcousticDriver {
    AcousticWorkspace workspace_;
    Spectrum span_, lengths_;
    std::size_t blades_, first_;
    std::vector<Vec3> observers_;
    double dt_, start_, last_time_ = -std::numeric_limits<double>::infinity();
    std::vector<Node> selected_;
    Spectrum speeds_;
    std::vector<Vec3> inflow_, leading_;

    TurbulenceState state_;
    const Snapshot *advance(double, const std::vector<std::vector<Node>> &, const ObserverBlockCallback &,
                            std::size_t);

  public:
    const TurbulenceState &turbulence_state() const noexcept { return state_; }
    void set_propagation(std::shared_ptr<const OutdoorPropagation> p) {
        workspace_.set_propagation(std::move(p));
    }
    AcousticDriver(Parameters, Spectrum span, std::size_t blades, std::vector<Vec3> observers, double dt = .1,
                   double start = 0., double percentage = 70., double hub_height = 0., int ti_method = 1);
    bool is_sample_time(double time) const;
    std::size_t first_node() const noexcept { return first_; }
    std::optional<Snapshot> step(double time, const std::vector<std::vector<Node>> &blades);
    // Borrowed output; nullptr on non-sampling steps. Valid until the next call.
    const Snapshot *step_view(double time, const std::vector<std::vector<Node>> &blades);
    bool step_blocks(double time, const std::vector<std::vector<Node>> &blades, const ObserverBlockCallback &,
                     std::size_t block_size = 1);
};
std::vector<Vec3> read_observers(const std::string &path);
std::pair<Parameters, std::map<std::string, std::string>> read_aa_input(const std::string &path);
std::string backend();
void validate(const Parameters &);
struct QuadratureResult {
    double value, error, absolute, ascending;
};
QuadratureResult qk61(const std::function<double(double)> &, double lower, double upper);
double dot(const Spectrum &, const Spectrum &);
Spectrum a_weighting(const Spectrum &);
double db_sum(const Spectrum &);
// Legacy raw kernel samples. Arbitrary increasing frequencies remain accepted here;
// they do not establish non-overlapping bands or a PSD definition.
Mechanisms section_spectrum(const Parameters &, const Section &);
std::array<BandSoundPressureLevel, mechanism_count> band_section_spectrum(const Parameters &,
                                                                          const Section &);
std::pair<Spectrum, Spectrum> tblte_tno(double speed, double theta, double phi, double span, double distance,
                                        const BoundaryLayer &, const Parameters &);
double spl_integrate(double, double, double, bool, double, const BoundaryLayer &, const Parameters &);
std::pair<Geometry, Geometry> observe(const Vec3 &, const Vec3 &, const Mat3 &, double,
                                      std::array<double, 2> reference = {{.25, 0.}});
} // namespace aeroacoustics
