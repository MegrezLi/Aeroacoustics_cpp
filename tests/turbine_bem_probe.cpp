#include "turbine/bem.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
using namespace turbine;
int main(int argc,char** argv){try{
    if(argc!=4)throw std::runtime_error("Usage: turbine_bem_probe CASE.fst FIXTURE OUTPUT.csv");
    Case c(argv[1]);std::ofstream fixture(argv[2]),out(argv[3]);fixture<<std::setprecision(17);out<<std::setprecision(17)<<"sample,residual,a,ap,k,kp,F,valid\n";
    fixture<<c.airfoils.size()<<' '<<7200<<'\n';for(const auto& af:c.airfoils)fixture<<std::filesystem::absolute(af.input.path).generic_string()<<'\n';
    for(int n=0;n<7200;++n){const auto& s=c.stations[n%30];BEMOptions p;p.axial_drag=(n/30)%2;p.tangential_drag=(n/60)%2;p.tangential=(n/120)%2;
        BEMInput u{2+s.span,s.chord,s.twist,6+2*std::sin(.19*n),5+60*(.5+.5*std::sin(.13*n)),.4+std::abs(std::sin(.03*n)),.1+2*std::abs(std::sin(.07*n))};
        const double phi=(n/240)%2?-.7-.4*std::sin(.1*n):.35+.3*std::sin(.1*n);if((n/480)%2)u.vx=-u.vx;
        const auto r=bem_residual(p,u,c.airfoils[s.airfoil],phi);
        fixture<<s.airfoil+1<<' '<<p.axial_drag<<' '<<p.tangential_drag<<' '<<p.tangential<<' '<<u.radius<<' '<<u.chord<<' '<<u.twist<<' '<<u.vx<<' '<<u.vy<<' '<<u.hub_loss_constant<<' '<<u.tip_loss_constant<<' '<<phi<<'\n';
        out<<n<<','<<r.residual<<','<<r.axial<<','<<r.tangential<<','<<r.k<<','<<r.kp<<','<<r.loss<<','<<r.valid<<'\n';}
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
