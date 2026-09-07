// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include "kernels.hpp"
#include "quadrature_data.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#ifdef AERO_USE_MKL
#include <mkl.h>
#endif

namespace aeroacoustics {
namespace {
constexpr double pi=3.14159265358979323846;
void exponential(const Spectrum& x,Spectrum& y) {
    y.resize(x.size());
#ifdef AERO_USE_MKL
    vdExp(static_cast<MKL_INT>(x.size()),x.data(),y.data());
#else
    std::transform(x.begin(),x.end(),y.begin(),[](double a){return std::exp(a);});
#endif
}
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
} // namespace aeroacoustics
