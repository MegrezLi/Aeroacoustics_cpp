// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include "kernels.hpp"
#include "quadrature_data.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#ifdef AERO_USE_MKL
#include <mkl.h>
#endif

namespace aeroacoustics {
namespace {
constexpr double pi=3.14159265358979323846;
constexpr double silence=-std::numeric_limits<double>::infinity();
void exponential(const Spectrum& x,Spectrum& y) {
    y.resize(x.size());
#ifdef AERO_USE_MKL
    vdExp(static_cast<MKL_INT>(x.size()),x.data(),y.data());
#else
    std::transform(x.begin(),x.end(),y.begin(),[](double a){return std::exp(a);});
#endif
}
}
std::string backend() {
#ifdef AERO_USE_MKL
    char text[256]{};mkl_get_version_string(text,256);
    return std::string("Intel oneMKL: ")+text;
#else
    return "Portable C++17";
#endif
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
QuadratureResult qk61(const std::function<double(double)>& function,double lower,double upper) {
    const double centre=(lower+upper)/2,half=(upper-lower)/2,fc=function(centre);
    double resk=wgk[30]*fc,resg=0,resabs=wgk[30]*std::abs(fc);
    std::array<double,30> fm{},fp{};
    for(int i=0;i<30;++i) {
        fm[i]=function(centre-half*xgk[i]);fp[i]=function(centre+half*xgk[i]);
        resk+=wgk[i]*(fm[i]+fp[i]);resabs+=wgk[i]*(std::abs(fm[i])+std::abs(fp[i]));
        if(i%2==1)resg+=wg[i/2]*(fm[i]+fp[i]);
    }
    const double mean=resk/2;double resasc=wgk[30]*std::abs(fc-mean);
    for(int i=0;i<30;++i)resasc+=wgk[i]*(std::abs(fm[i]-mean)+std::abs(fp[i]-mean));
    resabs*=std::abs(half);resasc*=std::abs(half);double error=std::abs((resk-resg)*half);
    if(resasc!=0&&error!=0)error=resasc*std::min(1.,std::pow(200*error/resasc,1.5));
    const double eps=std::numeric_limits<double>::epsilon(),tiny=std::numeric_limits<double>::min();
    if(resabs>tiny/(50*eps))error=std::max(50*eps*resabs,error);
    return {resk*half,error,resabs,resasc};
}
double dot(const Spectrum& x,const Spectrum& y) {
    if(x.size()!=y.size())throw std::invalid_argument("dot: length mismatch");
#ifdef AERO_USE_MKL
    return cblas_ddot(static_cast<MKL_INT>(x.size()),x.data(),1,y.data(),1);
#else
    return std::inner_product(x.begin(),x.end(),y.begin(),0.);
#endif
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
double spl_integrate(double omega,double lower,double upper,bool suction,double mach,
                     const BoundaryLayer& bl,const Parameters& p) {
    const int side=suction?0:1;
    const double cf=bl.cf[side],delta=bl.d99[side];
    const double uo=mach*p.spdsound*std::abs(bl.edge_velocity_ratio[side]);
    if(cf<=0||delta<=0||uo<=0)throw std::invalid_argument("TNO needs positive Cf, delta and edge speed");
    Spectrum nodes(61),weights(61),wave(61),x2(61),factor(61*61),gauss(61*61),decay(61*61);
    for(int i=0;i<30;++i){nodes[i]=-xgk[i];nodes[i+31]=xgk[i];weights[i]=weights[i+31]=wgk[i];}
    nodes[30]=0.;weights[30]=wgk[30];
    for(int i=0;i<61;++i){wave[i]=(lower+upper)/2+(upper-lower)/2*nodes[i];x2[i]=delta/2*(1+nodes[i]);}
    const double alpha=suction?.45:.30,kappa=.41,cmu=.09,cnuk=5.5;
    const double ustar=uo*std::sqrt(cf/2.);
    const double wf=uo/ustar-std::log(ustar*delta/p.kinvisc)/kappa-cnuk;
    for(int j=0;j<61;++j) {
        const double x=x2[j],length=.085*delta*std::tanh(kappa*x/(.085*delta));
        const double u=ustar*(std::log(ustar*x/p.kinvisc)/kappa+cnuk+wf*.5*(1-std::cos(pi*x/delta)));
        const double grad=ustar*(1/(kappa*x)+wf*.5*pi/delta*std::sin(pi*x/delta));
        const double ke=std::sqrt(pi)/length*.4213560764;
        const double nut=std::pow(length*kappa,2)*std::abs(grad);
        const double ums=alpha*std::sqrt(std::pow(nut*grad,2)/cmu),uc=.7*u,ag=.05*uc/length;
        for(int i=0;i<61;++i) {
            const int idx=i*61+j;const double kh=wave[i]/ke;
            const double phi22=4./9./pi/(ke*ke)*(kh*kh)/std::pow(1+kh*kh,7./3.);
            factor[idx]=length*ums*grad*grad*phi22/(ag*std::sqrt(pi));
            gauss[idx]=-std::pow((omega-uc*wave[i])/ag,2);
            decay[idx]=-2*std::abs(wave[i])*x;
        }
    }
    Spectrum e1,e2;exponential(gauss,e1);exponential(decay,e2);
    for(std::size_t i=0;i<factor.size();++i)factor[i]*=e1[i]*e2[i];
    Spectrum pressure(61,0.);
#ifdef AERO_USE_MKL
    cblas_dgemv(CblasRowMajor,CblasNoTrans,61,61,1.,factor.data(),61,weights.data(),1,0.,pressure.data(),1);
#else
    for(int i=0;i<61;++i)for(int j=0;j<61;++j)pressure[i]+=factor[i*61+j]*weights[j];
#endif
    for(int i=0;i<61;++i)pressure[i]*=4*p.airdens*p.airdens*delta/2*omega/p.spdsound/wave[i];
    return dot(weights,pressure)*(upper-lower)/2;
}
std::pair<Spectrum,Spectrum> tblte_tno(double u,double theta,double phi,double span,double r,
                                     const BoundaryLayer& bl,const Parameters& p) {
    Spectrum suction(p.freqlist.size(),0.),pressure=suction;
    const double mach=u/p.spdsound,directivity=directh_te(mach,theta,phi),ratio=std::pow(2.,1./3.);
    for(std::size_t i=0;i<p.freqlist.size();++i) {
        const double omega=2*pi*p.freqlist[i];
        const double width=2*omega*(std::sqrt(ratio)-1/std::sqrt(ratio));
        for(int side=0;side<2;++side) {
            double& value=side==0?suction[i]:pressure[i];
            if(bl.cf[side]>0) {
                const double integral=spl_integrate(omega,0,10*omega/u,side==0,mach,bl,p);
                value=10*std::log10(span/(4*pi*r*r)*integral*directivity/(2e-5*2e-5))+10*std::log10(width);
            }
            value=std::max(value,-100.);
        }
    }
    return {pressure,suction};
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
std::pair<Geometry,Geometry> observe(const Vec3& obs,const Vec3& center,const Mat3& rot,double chord,std::array<double,2> ref) {
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){double value=0;
        for(int k=0;k<3;++k)value+=rot[i*3+k]*rot[j*3+k];
        if(!std::isfinite(value)||std::abs(value-(i==j?1.:0.))>1e-8+(i==j?1e-5:0.))
            throw std::invalid_argument("global_to_local must be orthonormal");}
    std::array<Geometry,2> output;
    for(int edge=0;edge<2;++edge) {
        Vec3 offset{{-ref[1]*chord,(edge-ref[0])*chord,0}},global{},local{};
        for(int i=0;i<3;++i){global[i]=obs[i]-center[i];for(int j=0;j<3;++j)global[i]-=offset[j]*rot[j*3+i];}
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)local[i]+=rot[i*3+j]*global[j];
        const double phi=std::atan2(local[0],local[2]);
        const double theta=std::atan2(local[2]*std::cos(phi)+local[0]*std::sin(phi),local[1]);
        output[edge]={std::max(1e-16,std::sqrt(local[0]*local[0]+local[1]*local[1]+local[2]*local[2])),theta*180/pi,phi*180/pi};
    }
    return {output[0],output[1]};
}
}
