// Copyright (C) 2012-2016 National Renewable Energy Laboratory.
// C++ derivative of OpenFAST, Apache-2.0; see LICENSE and NOTICE.
#include "aeroacoustics.hpp"
#include <numeric>
#include <stdexcept>
#ifdef AERO_USE_MKL
#include <mkl.h>
#endif

namespace aeroacoustics {
std::string backend() {
#ifdef AERO_USE_MKL
    char text[256]{};mkl_get_version_string(text,256);
    return std::string("Intel oneMKL: ")+text;
#else
    return "Portable C++17";
#endif
}
double dot(const Spectrum& x,const Spectrum& y) {
    if(x.size()!=y.size())throw std::invalid_argument("dot: length mismatch");
#ifdef AERO_USE_MKL
    return cblas_ddot(static_cast<MKL_INT>(x.size()),x.data(),1,y.data(),1);
#else
    return std::inner_product(x.begin(),x.end(),y.begin(),0.);
#endif
}
} // namespace aeroacoustics
