#include "aeroacoustics.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace aeroacoustics;
template<class T>void emit(const T& values){for(double v:values)std::cout<<v<<' ';}
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::cout<<std::setprecision(17);
    auto table=BLTable::read(argv[1]);std::cout<<"bl ";
    for(double a:{-30.,-2.,15.,40.})for(double r:{.1e6,2e6,9e6}) {
        auto bl=table.interpolate(a,r,2.7);emit(bl.dstar);emit(bl.d99);emit(bl.cf);emit(bl.edge_velocity_ratio);
    }
    std::cout<<"\nelements ";
    for(double percent:{20.,70.,100.}){auto e=blade_elements({0,1,5,10,16},percent);std::cout<<e.first<<' ';emit(e.second);}
    std::cout<<"\nquadrature ";
    for(double upper:{-1.,0.,.5,3.}){auto q=qk61([](double x){return std::exp(x)*std::cos(3*x);},-.3,upper);emit(std::array<double,4>{q.value,q.error,q.absolute,q.ascending});}
    std::cout<<"\nthickness ";emit(guidati_thickness({{.25,0},{1,0},{.1,-.05},{.01,-.02},{0,0},{.01,.03},{.1,.07},{1,0}}));
    Parameters p;p.itrip=0;p.lammod=1;p.tipmod=1;p.bluntmod=1;p.timod=2;p.aweighting=true;
    for(int method:{1,2}) {
        AcousticDriver driver(p,{1,5,10},2,{{175,0,2},{0,175,2}},1.,0.,70.,10.,method);
        std::cout<<"\ndriver"<<method<<' ';
        for(int step=0;step<9;++step) {
            std::vector<std::vector<Node>> blades(2,std::vector<Node>(3));
            for(int b=0;b<2;++b)for(int j=0;j<3;++j) {
                auto& n=blades[b][j];n.section.speed=40+3*j+b;n.section.chord=.5+.2*j;
                n.section.alpha_deg=2+j;n.aero_center={{0.,double(b*2),10.+(j+1)*3}};
                n.inflow={{8.+std::sin(step+b+j),.2*step,0}};
                double a=.1*(step+b);n.global_to_local={{std::cos(a),-std::sin(a),0,std::sin(a),std::cos(a),0,0,0,1}};
            }
            auto result=driver.step(step,blades);
            if(!result)return 3;
            for(const auto& observer:*result)for(const auto& node:observer)for(const auto& mechanism:node)emit(mechanism);
            emit(driver.turbulence_state().values);
        }
    }
    std::cout<<'\n';
}
