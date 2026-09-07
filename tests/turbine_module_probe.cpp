#include "turbine/input.hpp"
#include "turbine/bem.hpp"
#include "turbine/unsteady.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        if(argc!=3)throw std::runtime_error("Usage: turbine_module_probe CASE.fst OUTPUT.csv");
        turbine::Case c(argv[1]);std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot create probe output");
        std::cout<<"Input parsed: "<<c.stations.size()<<" blade stations, "<<c.airfoils.size()<<" airfoils, "<<c.duration<<" seconds\n";
        out<<"foil,step,alpha_rad,speed,cl,cd,cm,phi,axial,tangential,residual\n"<<std::setprecision(17);
        turbine::BEMOptions options;
        const double omega=c.structure.number("RotSpeed")*2*turbine::pi/60,hub=c.structure.number("HubRad"),tip=c.structure.number("TipRad");
        for(std::size_t i=0;i<c.stations.size();++i) {
            const auto& s=c.stations[i];const auto& af=c.airfoils[s.airfoil];
            turbine::UnsteadyAirfoil ua(af,s.chord,c.dt,c.sound_speed);
            double phi=.1;
            for(std::size_t n=0;n<1200;++n) {
                const double t=n*c.dt,alpha=(8+6*std::sin(2*turbine::pi*.7*t))*turbine::deg,speed=45+5*std::sin(.8*t);
                const auto cf=ua.evaluate(alpha,speed);
                turbine::BEMInput in{hub+s.span,s.chord,s.twist+c.structure.number("BlPitch(1)")*turbine::deg,8,omega*(hub+s.span),3*s.span/(2*hub),3*(tip-hub-s.span)/(2*(hub+s.span))};
                const auto bem=turbine::solve_bem(options,in,af,phi);phi=bem.phi;
                if(!std::isfinite(cf.cl+cf.cd+cf.cm))throw std::runtime_error("Nonfinite UA output");
                out<<i<<','<<n<<','<<alpha<<','<<speed<<','<<cf.cl<<','<<cf.cd<<','<<cf.cm<<','<<phi<<','<<bem.axial<<','<<bem.tangential<<','<<bem.residual<<'\n';
                ua.advance(alpha,speed,n);
            }
        }
        std::cout<<"Module diagnostic completed; this is not the coupled turbine simulation.\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
