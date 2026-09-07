// C++ implementation of OpenFAST acoustic driver and auxiliary algorithms.
#include "aeroacoustics.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>
#include <numeric>

namespace aeroacoustics {
namespace {
constexpr double pi=3.14159265358979323846;
void require(bool condition,const char* message){if(!condition)throw std::invalid_argument(message);}
std::vector<std::string> lines(const std::string& path) {
    std::ifstream stream(path);require(bool(stream),"Cannot read input file");
    std::vector<std::string> output;std::string line;
    while(std::getline(stream,line))output.push_back(line);
    if(!output.empty()&&output[0].compare(0,3,"\xef\xbb\xbf")==0)output[0].erase(0,3);
    return output;
}
double norm(const Vec3& v){return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}
void increasing(const Spectrum& v) {
    require(!v.empty(),"Empty interpolation axis");
    for(std::size_t i=0;i<v.size();++i)
        require(std::isfinite(v[i])&&(i==0||v[i]>v[i-1]),"Axes must be finite and strictly increasing");
}
std::tuple<std::size_t,std::size_t,double> bounds(const Spectrum& grid,double x) {
    if(grid.size()==1)return {0,0,0};
    auto pos=std::upper_bound(grid.begin(),grid.end(),x)-grid.begin();
    std::size_t a=std::clamp<std::ptrdiff_t>(pos-1,0,grid.size()-2);
    return {a,a+1,std::clamp((x-grid[a])/(grid[a+1]-grid[a]),0.,1.)};
}
}
BLTable BLTable::read(const std::string& path) {
    auto text=lines(path);require(text.size()>=4,"Truncated BL table");
    int nr=std::stoi(text[2]),na=std::stoi(text[3]);require(nr>0&&na>0,"Invalid BL table dimensions");
    BLTable result;std::size_t pos=4;Spectrum raw_aoa;
    for(int r=0;r<nr;++r) {
        while(pos<text.size()&&(text[pos].find_first_not_of(" \t\r")==std::string::npos||
            text[pos][text[pos].find_first_not_of(" \t\r")]=='!'))++pos;
        require(pos+3+na<=text.size(),"Truncated BL rows");
        result.reynolds.push_back(std::stod(text[pos])*1e6);pos+=3;
        for(int a=0;a<na;++a) {
            std::istringstream row(text[pos++]);double angle;std::array<double,8> values{};
            require(bool(row>>angle)&&std::isfinite(angle),"Invalid BL AoA");
            for(auto& v:values)require(bool(row>>v)&&std::isfinite(v),"Invalid BL value");
            if(r==0){raw_aoa.push_back(angle);result.aoa.push_back(angle-360*std::floor((angle+180)/360));}
            else require(std::abs(raw_aoa[a]-angle)<=1e-8+1e-5*std::abs(raw_aoa[a]),"BL AoA grids differ");
            result.values.push_back(values);
        }
    }
    increasing(result.aoa);increasing(result.reynolds);return result;
}
BoundaryLayer BLTable::interpolate(double alpha,double re,double chord) const {
    increasing(aoa);increasing(reynolds);
    require(values.size()==aoa.size()*reynolds.size(),"Invalid BL shape");
    require(std::isfinite(alpha)&&std::isfinite(re)&&chord>0,"Invalid BL interpolation input");
    auto [a,b,t]=bounds(aoa,alpha);auto [c,d,s]=bounds(reynolds,re);
    std::array<double,8> v{};auto n=aoa.size();
    for(int k=0;k<8;++k)v[k]=(1-s)*((1-t)*values[c*n+a][k]+t*values[c*n+b][k])+
        s*((1-t)*values[d*n+a][k]+t*values[d*n+b][k]);
    return {{{v[2]*chord,v[3]*chord}},{{v[4]*chord,v[5]*chord}},{{v[6],v[7]}},{{v[0],v[1]}}};
}
std::pair<std::size_t,Spectrum> blade_elements(const Spectrum& span,double percent) {
    increasing(span);require(span.front()>=0&&percent>0&&percent<=100,"Invalid span or percentage");
    const auto n=span.size();std::size_t start=std::max<std::size_t>(1,n-1);
    double threshold=span.back()*(1-percent/100);
    for(std::size_t j=n-1;j>1;--j)if(span[j-1]<threshold){start=j;break;}
    start=std::max(std::min<std::size_t>(n,2),start)-1;
    Spectrum length(n,0);
    for(auto j=start;j<n;++j)length[j]=j==0?span[j]:j==n-1?span[j]-span[j-1]:(span[j+1]-span[j-1])/2;
    return {start,length};
}
std::array<double,2> guidati_thickness(const std::vector<std::array<double,2>>& coords) {
    require(coords.size()>=3,"Airfoil needs reference and contour rows");
    std::size_t split=coords.size();
    for(std::size_t j=2;j<coords.size();++j)if(coords[j][0]>coords[j-1][0]){split=j;break;}
    require(split<coords.size(),"Expected TE to LE to TE airfoil contour");
    std::array<double,2> output{};int k=0;
    for(double target:{.01,.10}) {
        std::size_t upper=0,lower=split;
        for(std::size_t j=1;j<split;++j)if(std::abs(coords[j][0]-target)<std::abs(coords[upper][0]-target))upper=j;
        for(std::size_t j=split+1;j<coords.size();++j)if(std::abs(coords[j][0]-target)<std::abs(coords[lower][0]-target))lower=j;
        output[k++]=coords[lower][1]-coords[upper][1];
    }
    return output;
}
TurbulenceState::TurbulenceState(Spectrum span,std::size_t blades,double dt,double height,int method,double ti,double avgv):
    span_(std::move(span)),blades_(blades),height_(height),ti_(ti),avgv_(avgv),method_(method) {
    increasing(span_);require(span_[0]>=0&&blades>0&&dt>0&&(method==1||method==2),"Invalid TI initialization");
    samples_=std::max<std::size_t>(std::floor(5/dt+.5),1);
    auto regions=static_cast<std::size_t>(std::ceil(span_.back()/5))+1;
    for(double r:span_)radial_.push_back(std::clamp<std::size_t>(std::ceil(r/5),1,regions)-1);
    buffers_.assign(regions*6*samples_,0);counts_.assign(regions*6,0);values.assign(blades*span_.size(),0);
}
void TurbulenceState::update(const Spectrum& vrel,const std::vector<Vec3>& inflow,const std::vector<Vec3>& leading) {
    require(vrel.size()==values.size()&&inflow.size()==values.size()&&leading.size()==values.size(),"TI array shape mismatch");
    for(std::size_t k=0;k<values.size();++k) {
        if(method_==1){require(vrel[k]!=0,"TI scaling requires nonzero Vrel");values[k]=ti_*avgv_/vrel[k];continue;}
        double z=leading[k][2]-height_;std::size_t angular=0;
        if(std::abs(z)>std::max(std::abs(z),1.)*std::numeric_limits<double>::epsilon()*50) {
            double angle=std::atan2(leading[k][1],z)*180/pi;angle-=360*std::floor(angle/360);
            angular=std::clamp<int>(std::ceil(angle/60),1,6)-1;
        }
        const auto region=radial_[k%span_.size()]*6+angular,count=++counts_[region];
        auto begin=buffers_.begin()+region*samples_;
        if(count<=samples_){begin[count-1]=norm(inflow[k]);values[k]=0;}
        else {
            begin[count%samples_]=norm(inflow[k]);
            double mean=std::accumulate(begin,begin+samples_,0.)/samples_,variance=0;
            for(std::size_t j=0;j<samples_;++j)variance+=std::pow(begin[j]-mean,2);
            values[k]=std::abs(mean)>std::max(std::abs(mean),1.)*std::numeric_limits<double>::epsilon()*50?
                std::sqrt(variance/samples_)/mean:0.;
        }
    }
}
Snapshot snapshot_spectrum(const Parameters& p,const std::vector<Node>& nodes,const std::vector<Vec3>& observers) {
    Snapshot result;
    for(const auto& observer:observers) {
        std::vector<Mechanisms> spectra;
        for(const auto& node:nodes) {
            auto s=node.section;auto geometry=observe(observer,node.aero_center,node.global_to_local,s.chord,node.airfoil_reference);
            s.leading=geometry.first;s.trailing=geometry.second;spectra.push_back(section_spectrum(p,s));
        }
        result.push_back(std::move(spectra));
    }
    return result;
}
AcousticDriver::AcousticDriver(Parameters p,Spectrum span,std::size_t blades,std::vector<Vec3> observers,
    double dt,double start,double percentage,double height,int method):parameters_(p),span_(span),blades_(blades),
    observers_(std::move(observers)),dt_(dt),start_(start),state(span,blades,dt,height,method,p.ti,p.avgv) {
    std::tie(first_,lengths_)=blade_elements(span,percentage);require(!observers_.empty(),"No observers");
}
std::optional<Snapshot> AcousticDriver::step(double time,const std::vector<std::vector<Node>>& blades) {
    require(std::isfinite(time)&&time>last_time_,"Times must be finite and increasing");
    require(blades.size()==blades_,"Blade count mismatch");
    std::vector<Node> selected;Spectrum vrel;std::vector<Vec3> inflow,leading;
    for(std::size_t b=0;b<blades_;++b) {
        require(blades[b].size()==span_.size(),"Node count mismatch");
        for(std::size_t j=0;j<span_.size();++j) {
            auto node=blades[b][j];auto& s=node.section;vrel.push_back(s.speed);inflow.push_back(node.inflow);
            Vec3 position=node.aero_center;
            for(int k=0;k<3;++k)position[k]+=s.chord*(-node.airfoil_reference[1]*node.global_to_local[k]
                -node.airfoil_reference[0]*node.global_to_local[3+k]);
            leading.push_back(position);
            if(j>=first_){s.span=lengths_[j];s.is_tip=j==span_.size()-1;s.ti_section=state.values[b*span_.size()+j];selected.push_back(node);}
        }
    }
    std::optional<Snapshot> result;
    double phase=time+1e-10-dt_*std::floor((time+1e-10)/dt_);
    if(time>=start_&&phase<1e-6)result=snapshot_spectrum(parameters_,selected,observers_);
    state.update(vrel,inflow,leading);last_time_=time;return result;
}
std::vector<Vec3> read_observers(const std::string& path) {
    auto text=lines(path);require(text.size()>=3,"Truncated observer file");int n=std::stoi(text[0]);
    require(n>0&&text.size()>=std::size_t(n+2),"Invalid observer count");std::vector<Vec3> result(n);
    for(int i=0;i<n;++i){std::istringstream row(text[i+2]);for(auto& v:result[i])require(bool(row>>v)&&std::isfinite(v),"Invalid observer position");}
    return result;
}
std::pair<Parameters,std::map<std::string,std::string>> read_aa_input(const std::string& path) {
    std::map<std::string,std::string> fields;
    std::regex pattern(R"aa(^\s*("[^"]*"|'[^']*'|\S+)\s+([A-Za-z]\w*)\b)aa");
    auto lower=[](std::string s){for(auto& c:s)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));return s;};
    for(const auto& line:lines(path)) {
        std::smatch match;if(!std::regex_search(line,match,pattern))continue;
        std::string value=match[1];if(value.size()>1&&(value[0]=='"'||value[0]=='\''))value=value.substr(1,value.size()-2);
        fields[lower(match[2])]=value;
    }
    Parameters p;
    for(auto [key,member]:std::vector<std::pair<std::string,double Parameters::*>>{
        {"spdsound",&Parameters::spdsound},{"kinvisc",&Parameters::kinvisc},{"airdens",&Parameters::airdens},
        {"lturb",&Parameters::lturb},{"alprat",&Parameters::alprat},{"ti",&Parameters::ti},{"avgv",&Parameters::avgv}})
        if(fields.count(key))p.*member=std::stod(fields[key]);
    for(auto [key,member]:std::vector<std::pair<std::string,int Parameters::*>>{
        {"blmod",&Parameters::x_blmethod},{"tripmod",&Parameters::itrip},{"timod",&Parameters::timod},
        {"tbltemod",&Parameters::tbltemod},{"lammod",&Parameters::lammod},{"tipmod",&Parameters::tipmod},{"bluntmod",&Parameters::bluntmod}})
        if(fields.count(key))p.*member=std::stoi(fields[key]);
    for(auto [key,member]:std::vector<std::pair<std::string,bool Parameters::*>>{
        {"roundedtip",&Parameters::round},{"aweighting",&Parameters::aweighting}})if(fields.count(key)) {
            auto v=lower(fields[key]);require(v=="true"||v=="false","Invalid boolean in AA input");p.*member=v=="true";
        }
    return {p,fields};
}
}
