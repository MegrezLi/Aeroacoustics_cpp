// Compile this same file against baseline and candidate headers/libraries.
#include "aeroacoustics.hpp"
#include "propagation.hpp"
#include "turbine/metrics.hpp"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <new>
static bool counting=false;
static std::size_t allocations=0, bytes=0;
void* operator new(std::size_t n) {
    void* p=std::malloc(n ? n : 1); if(!p) throw std::bad_alloc();
    if(counting) { ++allocations; bytes+=n; } return p;
}
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p,std::size_t) noexcept { std::free(p); }
void operator delete[](void* p,std::size_t) noexcept { std::free(p); }
using namespace aeroacoustics;
volatile double checksum=0;
template<class F> void measure(const char* label,int repeats,F f,bool warm=true) {
    if(warm) f();
    allocations=bytes=0; counting=true;
    const auto start=std::chrono::steady_clock::now();
    for(int i=0;i<repeats;++i) f();
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    counting=false;
    std::cout<<label<<','<<repeats<<','<<allocations<<','<<bytes<<','<<seconds<<','<<checksum<<'\n';
}
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::invalid_argument("Usage: optimization_benchmark OUTPUT_DIR");
        const std::filesystem::path directory=argv[1];
        std::filesystem::create_directories(directory);
        Parameters p; p.timod=0; p.itrip=0; p.lammod=1;
        std::vector<Node> input(1); input[0].aero_center={0,0,80};
        std::vector<Vec3> receivers{{100,20,2}};
        auto s=input[0].section;
        auto g=observe(receivers[0],input[0].aero_center,input[0].global_to_local,s.chord,input[0].airfoil_reference);
        s.leading=g.first; s.trailing=g.second;
        AcousticWorkspace workspace(p);
        std::cout<<std::setprecision(17)<<"case,repeats,allocations,requested_bytes,seconds,checksum\n";
        measure("BPM_section",1000,[&]{auto result=section_spectrum(p,s);checksum+=result[0][10];});
        measure("BPM_workspace",1000,[&]{const auto& result=workspace.evaluate(input,receivers);checksum+=result[0][0][0][10];});
        input.resize(66);
        for(std::size_t i=0;i<input.size();++i) input[i].aero_center={0,double(i),80};
        receivers.assign(17,{100,20,2});
        AcousticWorkspace block(p);
        ObserverBlockCallback callback=[](std::size_t,const Snapshot& result){checksum+=result[0][0][0][10];};
        measure("17_observers_block4",30,[&]{block.evaluate_blocks(input,receivers,callback,4);});

        PropagationOptions propagation;
        propagation.ground=GroundModel::impedance;
        propagation.screens={{20,-100,20,100,50},{50,-100,50,100,30}};
        OutdoorPropagation outdoor(propagation,FrequencyBands::openfast_reference(p.freqlist));
        Mechanisms direct, image, working;
        for(std::size_t m=0;m<mechanism_count;++m) {
            direct[m].assign(p.freqlist.size(),60.+m);
            image[m].assign(p.freqlist.size(),55.+m);
        }
        working=direct;
        measure("propagation",10000,[&]{working=direct;outdoor.apply(input[0],receivers[0],working,&image);checksum+=working[0][10];});

        TurbulenceState ti({5.},1,.005,80.,2);
        Spectrum speed{50.};
        std::vector<Vec3> wind{{8,0,0}},position{{0,0,80}};
        for(int i=0;i<1100;++i) ti.update(speed,wind,position);
        int iteration=0;
        measure("TI_method2",100000,[&]{wind[0][0]=8.+std::sin(++iteration*.01);ti.update(speed,wind,position);checksum+=ti.values[0];});

        const std::size_t nodes=12, observers=8;
        receivers.resize(observers);
        for(std::size_t o=0;o<observers;++o) receivers[o]={100.+o,20.,2.};
        std::vector<std::vector<Node>> blades(1,std::vector<Node>(nodes));
        for(std::size_t n=0;n<nodes;++n) blades[0][n].aero_center={0.,double(n),80.};
        Snapshot snapshot(observers,std::vector<Mechanisms>(nodes,direct));
        turbine::MetricsOptions options;
        options.start=1.; options.end=19.; options.receiver_dt=.1;
        turbine::EngineeringMetrics metrics(options,p,receivers,nodes,{0,0,80},std::nullopt);
        int sample=0;
        measure("metrics_append",201,[&]{metrics.append(.1*sample++,8.,1.,blades,0,0,snapshot);},false);
        measure("metrics_write",1,[&]{metrics.write(directory);},false);
    } catch(const std::exception &e) {
        counting=false;
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
