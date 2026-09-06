#pragma once
#include <array>
#include <string>
#include <tuple>
#include <vector>
#include <optional>
#include <limits>
#include <map>
#include <functional>

namespace aeroacoustics {
using Spectrum = std::vector<double>;
using Mechanisms = std::array<Spectrum,7>;
using Vec3 = std::array<double,3>;
using Mat3 = std::array<double,9>; // row-major global-to-local
struct Parameters {
    Spectrum freqlist = {10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000};
    double spdsound=340.,kinvisc=1.48e-5,airdens=1.225,lturb=40.,alprat=1.,ti=.1,avgv=8.;
    int x_blmethod=1,itrip=1,timod=1,tbltemod=1,lammod=0,tipmod=0,bluntmod=0;
    bool round=true,aweighting=false;
};
struct BoundaryLayer {
    std::array<double,2> dstar{},d99{},cf{},edge_velocity_ratio{{1.,1.}};
};
struct Geometry {double distance=1.22,theta=90.,phi=90.;};
struct Section {
    double chord=.2286,speed=63.92,alpha_deg=3.,span=.509,stall_deg=12.5;
    Geometry trailing{};
    std::optional<Geometry> leading;
    BoundaryLayer bl{};
    double ti_section=-1.,te_thickness=.001,te_angle=14.,thickness_1p=.02,thickness_10p=.12;
    bool is_tip=true;
};
struct BLTable {
    Spectrum aoa,reynolds;
    std::vector<std::array<double,8>> values; // Re-major, then AoA
    static BLTable read(const std::string& path);
    BoundaryLayer interpolate(double alpha_deg,double re,double chord) const;
};
struct Node {
    Section section;
    Vec3 aero_center{},inflow{};
    Mat3 global_to_local{{1,0,0,0,1,0,0,0,1}};
    std::array<double,2> airfoil_reference{{.25,0}};
};
using Snapshot=std::vector<std::vector<Mechanisms>>; // observer, node, mechanism, frequency
std::pair<std::size_t,Spectrum> blade_elements(const Spectrum& span,double percentage=100.);
std::array<double,2> guidati_thickness(const std::vector<std::array<double,2>>& coords);
Snapshot snapshot_spectrum(const Parameters&,const std::vector<Node>&,const std::vector<Vec3>&);
class TurbulenceState {
    Spectrum span_,buffers_;
    std::vector<std::size_t> radial_,counts_;
    std::size_t samples_,blades_;
    double height_,ti_,avgv_;
    int method_;
public:
    Spectrum values; // blade-major, then radial node
    TurbulenceState(Spectrum span,std::size_t blades,double dt,double hub_height,
                    int method=1,double ti=.1,double avgv=8.);
    void update(const Spectrum& vrel,const std::vector<Vec3>& inflow,const std::vector<Vec3>& leading);
};
class AcousticDriver {
    Parameters parameters_;
    Spectrum span_,lengths_;
    std::size_t blades_,first_;
    std::vector<Vec3> observers_;
    double dt_,start_,last_time_=-std::numeric_limits<double>::infinity();
public:
    TurbulenceState state;
    AcousticDriver(Parameters, Spectrum span,std::size_t blades,std::vector<Vec3> observers,
                   double dt=.1,double start=0.,double percentage=70.,double hub_height=0.,int ti_method=1);
    std::optional<Snapshot> step(double time,const std::vector<std::vector<Node>>& blades);
};
std::vector<Vec3> read_observers(const std::string& path);
std::pair<Parameters,std::map<std::string,std::string>> read_aa_input(const std::string& path);
std::string backend();
void validate(const Parameters&);
struct QuadratureResult {double value,error,absolute,ascending;};
QuadratureResult qk61(const std::function<double(double)>&,double lower,double upper);
double dot(const Spectrum&,const Spectrum&);
Spectrum a_weighting(const Spectrum&);
double db_sum(const Spectrum&);
Mechanisms section_spectrum(const Parameters&,const Section&);
std::pair<Spectrum,Spectrum> tblte_tno(double,double,double,double,double,const BoundaryLayer&,const Parameters&);
double spl_integrate(double,double,double,bool,double,const BoundaryLayer&,const Parameters&);
std::pair<Geometry,Geometry> observe(const Vec3&,const Vec3&,const Mat3&,double,std::array<double,2> reference = {{.25,0.}});
}
