#include "acoustic_levels.hpp"
#include "aeroacoustics_c.h"
#include "turbine/metrics.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace aeroacoustics;
namespace {
void check(bool ok, const char *message) {
    if (!ok) throw std::runtime_error(message);
}
void near(double a, double b) {
    check(std::isfinite(a) && std::isfinite(b) &&
          std::abs(a-b) <= 1e-10 * std::max({1., std::abs(a), std::abs(b)}), "Numeric mismatch");
}
template<class F> void rejects(F f) {
    bool failed=false;
    try { f(); } catch(const std::exception &) { failed=true; }
    check(failed, "Expected rejection");
}
void sums() {
    for (auto invalid : {NAN, INFINITY})
        for (Spectrum s : {Spectrum{-INFINITY, invalid}, Spectrum{invalid,-INFINITY},
                           Spectrum{20,invalid,30}})
            rejects([&] { db_sum(s); });
    check(db_sum({}) == -INFINITY && db_sum({-INFINITY,-INFINITY}) == -INFINITY, "Silence sum");
    near(db_sum({0,0,-INFINITY}),10*std::log10(2.));
    near(db_sum({4000,4000}),4000+10*std::log10(2.));
}
void percentiles() {
    for(double peak : {1e-50,100.,1e100}) {
        auto s=turbine::level_statistics({0,9,10},{0,0,peak});
        check(s.l50_db == -INFINITY && s.l95_db == -INFINITY, "Partial silence percentiles");
        near(s.l5_db,10*std::log10(.5*peak));
        near(s.leq_db,10*std::log10(.05*peak));
    }
    auto boundary=turbine::level_statistics({0,.5,10},{0,0,100});
    check(boundary.l95_db == -INFINITY, "Silence mass exactly at quantile");
    auto below=turbine::level_statistics({0,.4,10},{0,0,100});
    near(below.l95_db,10*std::log10(100.*.1/9.6));
    auto constant=turbine::level_statistics({0,1,2},{4,4,4});
    near(constant.l95_db,10*std::log10(4.));
    auto plateau=turbine::level_statistics({0,1,9,10},{0,4,4,8});
    near(plateau.l50_db,10*std::log10(4.));
    auto silence=turbine::level_statistics({0,1},{0,0});
    check(silence.leq_db == -INFINITY && silence.l5_db == -INFINITY &&
          silence.l50_db == -INFINITY && silence.l95_db == -INFINITY, "All-silent statistics");
}
void abi() {
    const std::array<double,26> valid{{3,.2,40,90,90,.5,10,12,.007,.001,.001,.1,.001,10,
                                      .02,.12,340,1.48e-5,1.225,40,1,.002,.002,.007,1,1}};
    const std::array<int,3> flags{{2,1,1}};
    const std::array<double,2> frequencies{{100,1000}};
    auto run=[&](int op, std::array<double,26> x, std::array<int,3> modes,
                 std::array<double,2> freq, bool fail) {
        std::array<double,6> out;
        out.fill(1234567.);
        auto before=aeroacoustics_call_count();
        int status=aeroacoustics_kernel(op,2,freq.data(),x.data(),modes.data(),out.data());
        check((status != 0) == fail, "Wrong C ABI status");
        if (fail) {
            check(std::string(aeroacoustics_last_error()).size()>0, "Missing C ABI diagnostic");
            check(aeroacoustics_call_count()==before, "Failed call incremented counter");
            for(double v:out) check(v==1234567.,"Failure partially wrote output");
        } else {
            check(std::string(aeroacoustics_last_error()).empty(), "Success retained stale error");
            check(aeroacoustics_call_count()==before+1, "Successful call counter");
            int channels=op==2?3:op==7?2:1;
            for(int i=0;i<2*channels;++i) validate_level(out[i]);
        }
        return out;
    };
    for(int op=1;op<=7;++op) {
        for(double v : {double(NAN),double(INFINITY),-1.,0.}) {
            auto f=frequencies; f[1]=v;
            run(op,valid,flags,f,true);
        }
        for(double v : {double(NAN),double(INFINITY),-1.,0.}) {
            auto x=valid; x[2]=v; run(op,x,flags,frequencies,true);
        }
        if(op!=6) {
            auto x=valid; x[6]=0; run(op,x,flags,frequencies,true);
            x=valid; x[2]=x[16]; run(op,x,flags,frequencies,true);
        }
        run(op,valid,flags,frequencies,false);
    }
    for(auto pair:std::vector<std::pair<int,int>>{{1,8},{2,9},{3,20},{4,11},{5,12},{6,14},{7,24}}) {
        auto x=valid; x[pair.second]=NAN;
        run(pair.first,x,flags,frequencies,true);
    }
    for(int op : {1,2,5}) {
        auto invalid=flags; invalid[0]=0;
        run(op,valid,invalid,frequencies,true);
    }
    auto invalid=flags; invalid[2]=2; run(3,valid,invalid,frequencies,true);
    auto x=valid; x[13]=15; run(5,x,flags,frequencies,true);
    x=valid; x[11]=-1; run(4,x,flags,frequencies,true);
    x=valid; x[11]=0; run(4,x,flags,frequencies,false);
    x=valid; x[20]=1e308; run(3,x,flags,frequencies,true); // finite inputs, invalid computed result
    x=valid;
    for(int i : {0,1,7,9,10,11,12,13,14,15,19,20}) x[i]=NAN;
    auto reference=run(7,x,flags,frequencies,false); // unused fields must not be inspected
    x[24]=-1; x[25]=-1;
    check(reference==run(7,x,flags,frequencies,false), "TNO signed ratio compatibility");
    x[24]=1.15; x[25]=.85;
    check(reference!=run(7,x,flags,frequencies,false), "TNO input ratios ignored");
    x[21]=0; x[22]=-1; x[24]=x[25]=0;
    x[8]=x[23]=0;
    run(7,x,flags,frequencies,false); // legacy nonpositive Cf disables a side
    x=valid; x[7]=x[9]=x[10]=NAN;
    run(1,x,{{2,-1,-1}},frequencies,false); // tabulated laminar kernel uses only d99P
    x.fill(NAN); x[1]=.2; x[2]=40; x[14]=.02; x[15]=.12;
    run(6,x,{{-1,-1,-1}},frequencies,false);
    std::array<double,6> out{};
    check(aeroacoustics_kernel(4,0,frequencies.data(),valid.data(),flags.data(),out.data())!=0,"Zero n");
    check(aeroacoustics_kernel(4,2,nullptr,valid.data(),flags.data(),out.data())!=0,"Null frequency");
    run(0,valid,flags,frequencies,true);
    run(4,valid,flags,frequencies,false);
}
std::vector<std::vector<double>> rows(const std::filesystem::path &file) {
    std::ifstream in(file);
    check(bool(in),"Missing CSV");
    std::string line,word;
    std::getline(in,line);
    std::vector<std::vector<double>> result;
    while(std::getline(in,line)) {
        std::istringstream stream(line);
        std::vector<double> row;
        while(std::getline(stream,word,',')) if(!word.empty()) row.push_back(std::stod(word));
        result.push_back(row);
    }
    return result;
}
double integral(const std::vector<double> &t,const std::vector<double> &p,double begin,double end) {
    double result=0;
    for(std::size_t i=1;i<t.size();++i) {
        double a=std::max(begin,t[i-1]), b=std::min(end,t[i]);
        if(a>=b) continue;
        auto at=[&](double x) { return p[i-1]+(p[i]-p[i-1])*(x-t[i-1])/(t[i]-t[i-1]); };
        result+=(b-a)*(at(a)+at(b))/2;
    }
    return result;
}
void metrics(const std::filesystem::path &output) {
    for(bool moving : {false,true}) {
        std::vector<double> reference;
        std::vector<std::vector<double>> reference_bins;
        for(double dt : {.1,.2,.07,.125}) {
            turbine::MetricsOptions options;
            options.start=moving?.12:0;
            options.end=moving?.36:.2;
            options.receiver_dt=dt;
            options.am_window=14;
            options.retarded_time=moving;
            options.background_laeq_db=0;
            Parameters p; p.freqlist={1000}; p.aweighting=true;
            const std::size_t n=moving?2:1;
            turbine::EngineeringMetrics m(options,p,{{100,0,0}},n,{100,0,0},{});
            std::vector<std::vector<Node>> nodes(1,std::vector<Node>(n));
            Snapshot snapshot(1,std::vector<Mechanisms>(n));
            for(auto &node:snapshot[0]) for(auto &source:node) source={-INFINITY};
            std::vector<std::vector<double>> times(n), powers(n);
            for(int i=0;i<(moving?5:3);++i) {
                double t=i*.1;
                for(std::size_t j=0;j<n;++j) {
                    double energy=j==0?(i==1?100:1):(i==2?50:2);
                    double delay=moving?(j==0?.05+.1*t:.1+.2*t):0;
                    nodes[0][j].airfoil_reference={1,0};
                    nodes[0][j].aero_center={100-p.spdsound*delay,0,0};
                    snapshot[0][j][1][0]=10*std::log10(energy);
                    times[j].push_back(t+delay); powers[j].push_back(energy);
                }
                m.append(t,i%2?8.5:7.5,1,nodes,0,0,snapshot);
            }
            auto folder=output/((moving?"moving-":"static-")+std::to_string(dt));
            std::filesystem::create_directories(folder);
            m.write(folder);
            auto summary=rows(folder/"receiver_map.csv").at(0);
            double energy=0;
            for(std::size_t j=0;j<n;++j)
                energy+=integral(times[j],powers[j],options.start,options.end);
            const double duration=options.end-options.start;
            near(summary[7],10*std::log10(energy/duration));
            near(summary[8],10*std::log10(energy/duration+1));
            auto bins=rows(folder/"wind_bins.csv");
            double bin_time=0,bin_energy=0;
            for(const auto &bin:bins) {
                bin_time+=bin[3]; bin_energy+=bin[3]*std::pow(10.,bin[4]/10);
            }
            near(bin_time,duration); near(bin_energy,energy);
            if(reference.empty()) { reference=summary; reference_bins=bins; }
            else {
                for(std::size_t k=4;k<summary.size();++k) near(summary[k],reference[k]);
                check(bins.size()==reference_bins.size(),"Wind bin count changed with ReceiverDT");
                for(std::size_t j=0;j<bins.size();++j)
                    for(std::size_t k=0;k<bins[j].size();++k) near(bins[j][k],reference_bins[j][k]);
            }
        }
    }
}
}
int main(int argc,char **argv) {
    try {
        check(argc==2,"review_regression_probe OUTPUT_DIRECTORY");
        sums(); percentiles(); abi(); metrics(argv[1]);
        std::cout<<"Four review regressions passed: energy conservation, C ABI validation, db_sum, silence quantiles\n";
    } catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; }
}
