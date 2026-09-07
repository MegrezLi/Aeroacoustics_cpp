// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace aeroacoustics {
namespace {
constexpr double pi=3.14159265358979323846;
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
} // namespace aeroacoustics
