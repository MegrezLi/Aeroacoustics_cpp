// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include "kernels.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aeroacoustics {
namespace {
constexpr double pi=3.14159265358979323846;
constexpr double silence=-std::numeric_limits<double>::infinity();
}
void validate(const Parameters& p) {
    if(p.freqlist.empty())throw std::invalid_argument("Empty frequency list");
    double last=0;
    for(double f:p.freqlist){if(!std::isfinite(f)||f<=last)throw std::invalid_argument("Frequencies must be finite, positive and increasing");last=f;}
    for(double v:{p.spdsound,p.kinvisc,p.airdens,p.lturb,p.alprat})
        if(!std::isfinite(v)||v<=0)throw std::invalid_argument("Positive physical parameters required");
    for(double v:{p.ti,p.avgv})if(!std::isfinite(v)||v<0)throw std::invalid_argument("Nonnegative TI and mean speed required");
    if(p.timod<0||p.timod>2||p.tbltemod<0||p.tbltemod>2||p.x_blmethod<1||p.x_blmethod>2||
       p.itrip<0||p.itrip>2||p.lammod<0||p.lammod>1||p.tipmod<0||p.tipmod>1||p.bluntmod<0||p.bluntmod>1)
        throw std::invalid_argument("Unsupported acoustic switch");
}
Spectrum a_weighting(const Spectrum& frequencies) {
    Spectrum result;result.reserve(frequencies.size());
    for(double f:frequencies) {
        const double f2=f*f,f4=f2*f2;
        result.push_back(10*std::log10(1.562339*f4/((f2+std::pow(107.65265,2))*(f2+std::pow(737.86223,2))))
           +10*std::log10(2.242881e16*f4/(std::pow(f2+std::pow(20.598997,2),2)*std::pow(f2+std::pow(12194.22,2),2))));
    }
    return result;
}
double db_sum(const Spectrum& levels) {
    if(levels.empty())return silence;
    const double maximum=*std::max_element(levels.begin(),levels.end());
    if(maximum==silence)return silence;
    double total=0.;for(double level:levels)total+=std::pow(10.,(level-maximum)/10.);
    return maximum+10*std::log10(total);
}
Mechanisms section_spectrum(const Parameters& input,const Section& s) {
    Parameters p=input;
    validate(p);
    for(double v:{s.speed,s.chord,s.span,s.alpha_deg,s.stall_deg,s.trailing.distance,s.trailing.theta,s.trailing.phi})
        if(!std::isfinite(v))throw std::invalid_argument("Section inputs must be finite");
    if(s.speed<0||s.speed>=p.spdsound||s.chord<=0||s.span<=0||s.trailing.distance<=0)
        throw std::invalid_argument("Invalid section speed or geometry");
    if(p.tbltemod==2)p.x_blmethod=2;
    if(p.x_blmethod==2&&(s.bl.dstar[0]<=0||s.bl.dstar[1]<=0||s.bl.d99[0]<=0||s.bl.d99[1]<=0))
        throw std::invalid_argument("Tabulated BL/TNO needs dimensional boundary-layer data");
    const double u=std::max(s.speed,.1),a=s.alpha_deg-360*std::floor((s.alpha_deg+180)/360);
    const auto& t=s.trailing;const auto le=s.leading.value_or(s.trailing);const auto& bl=s.bl;
    Mechanisms result;for(auto& v:result)v.assign(p.freqlist.size(),silence);
    if(p.lammod&&p.itrip==0)result[0]=lblvs(a,s.chord,u,t.theta,t.phi,s.span,t.distance,p,bl.d99[1],bl.dstar[0],bl.dstar[1],s.stall_deg);
    if(p.tbltemod) {
        std::tie(result[1],result[2],result[3])=tblte(a,s.chord,u,t.theta,t.phi,s.span,t.distance,p,bl.d99[1],bl.dstar[0],bl.dstar[1],s.stall_deg);
        if(p.tbltemod==2){auto one=bl;one.edge_velocity_ratio={{1,1}};auto v=tblte_tno(u,t.theta,t.phi,s.span,t.distance,one,p);result[1]=v.first;result[2]=v.second;}
    }
    if(p.bluntmod) {
        if(s.te_thickness<=0||s.te_angle<0||s.te_angle>14)throw std::invalid_argument("Invalid bluntness geometry");
        result[4]=blunt(a,s.chord,u,t.theta,t.phi,s.span,t.distance,s.te_thickness,s.te_angle,p,bl.d99[1],bl.dstar[0],bl.dstar[1],s.stall_deg);
    }
    if(p.tipmod&&s.is_tip)result[5]=tipnois(a,p.alprat,s.chord,u,t.theta,t.phi,t.distance,p);
    if(p.timod) {
        const double ti=s.ti_section<0?p.ti*p.avgv/u:s.ti_section;
        if(!std::isfinite(ti)||ti<0||le.distance<=0)throw std::invalid_argument("Invalid inflow noise input");
        result[6]=inflownoise(a*pi/180,s.chord,u,le.theta,le.phi,s.span,le.distance,ti,p);
        if(p.timod==2) {auto correction=simple_guidati(u,s.chord,s.thickness_10p,s.thickness_1p,p);
            for(std::size_t i=0;i<correction.size();++i)result[6][i]+=correction[i]+10.;}
        else if(p.timod!=1)throw std::invalid_argument("Unsupported TIMod");
    }
    if(p.aweighting){auto w=a_weighting(p.freqlist);for(auto& v:result)for(std::size_t i=0;i<v.size();++i)v[i]+=w[i];}
    return result;
}
} // namespace aeroacoustics
