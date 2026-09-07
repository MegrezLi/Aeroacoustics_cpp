// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include "quadrature_data.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace aeroacoustics {
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
} // namespace aeroacoustics
