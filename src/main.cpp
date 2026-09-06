#include "aeroacoustics.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    try {
        using namespace aeroacoustics;
        Parameters p;p.itrip=0;p.lammod=1;p.tipmod=1;p.bluntmod=1;
        Section section;
        auto result=section_spectrum(p,section);Spectrum total;
        const std::string filename=argc>1?argv[1]:"spectrum.csv";
        std::ofstream out(filename);if(!out)throw std::runtime_error("Cannot open output file");
        out<<"frequency_Hz,laminar,pressure,suction,separation,blunt,tip,inflow,total_dB\n"<<std::setprecision(16);
        for(std::size_t i=0;i<p.freqlist.size();++i){Spectrum components;out<<p.freqlist[i];
            for(const auto& v:result){out<<','<<v[i];components.push_back(v[i]);}
            total.push_back(db_sum(components));out<<','<<total.back()<<'\n';}
        std::cout<<backend()<<"\nOASPL = "<<std::setprecision(14)<<db_sum(total)<<" dB\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
