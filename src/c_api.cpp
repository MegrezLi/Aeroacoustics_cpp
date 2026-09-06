#include "aeroacoustics_c.h"
#include "aeroacoustics.hpp"
#include "kernels.hpp"
#include <algorithm>
#include <atomic>
#include <exception>
#include <stdexcept>
namespace {thread_local std::string error;std::atomic<unsigned long long> calls{0};}
extern "C" {
const char* aeroacoustics_last_error(){return error.c_str();}
const char* aeroacoustics_backend(){static std::string name=aeroacoustics::backend();return name.c_str();}
unsigned long long aeroacoustics_call_count(){return calls.load();}
int aeroacoustics_kernel(int op,int n,const double* f,const double* x,const int* flags,double* out) {
    try {
        using namespace aeroacoustics;
        if(n<=0||!f||!x||!flags||!out)throw std::invalid_argument("Invalid C ABI buffer");
        Parameters p;p.freqlist.assign(f,f+n);p.spdsound=x[16];p.kinvisc=x[17];p.airdens=x[18];p.lturb=x[19];
        p.x_blmethod=flags[0];p.itrip=flags[1];p.round=flags[2]!=0;
        Spectrum a,b,c;const auto copy=[&](const Spectrum& v,int k){std::copy(v.begin(),v.end(),out+k*n);};
        switch(op) {
        case 1:a=lblvs(x[0],x[1],x[2],x[3],x[4],x[5],x[6],p,x[8],x[9],x[10],x[7]);break;
        case 2:std::tie(a,b,c)=tblte(x[0],x[1],x[2],x[3],x[4],x[5],x[6],p,x[8],x[9],x[10],x[7]);copy(b,1);copy(c,2);break;
        case 3:a=tipnois(x[0],x[20],x[1],x[2],x[3],x[4],x[6],p);break;
        case 4:a=inflownoise(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[11],p);break;
        case 5:a=blunt(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[12],x[13],p,x[8],x[9],x[10],x[7]);break;
        case 6:a=simple_guidati(x[2],x[1],x[15],x[14],p);break;
        case 7:{BoundaryLayer bl;bl.cf={{x[21],x[22]}};bl.d99={{x[23],x[8]}};bl.edge_velocity_ratio={{x[24],x[25]}};
            auto pair=tblte_tno(x[2],x[3],x[4],x[5],x[6],bl,p);a=pair.first;copy(pair.second,1);break;}
        default:throw std::invalid_argument("Unknown acoustic kernel operation");
        }
        copy(a,0);++calls;error.clear();return 0;
    }catch(const std::exception& e){error=e.what();return 1;}catch(...){error="Unknown C++ failure";return 2;}
}
}
