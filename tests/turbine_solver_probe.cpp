#include "turbine/solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace turbine;
int main(int argc,char** argv){try{
    if(argc<3||argc>4)throw std::runtime_error("Usage: turbine_solver_probe CASE.fst OUTPUT.csv [duration]");
    Case c(argv[1]);if(argc==4)c.duration=std::stod(argv[3]);Solver solver(c);std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot create output");
    out<<std::setprecision(17)<<"time,blade,node,q1,q2,q3,qd1,qd2,qd3,phi,alpha,speed,a,ap,cl,cd,cm,x,y,z,fx,fy,fz,mx,my,mz\n";
    for(;;){for(int b=0;b<3;++b)for(std::size_t j=0;j<solver.aerodynamic().blades[b].size();++j){const auto& a=solver.aerodynamic().blades[b][j];out<<solver.time()<<','<<b<<','<<j;
        for(const auto& v:{solver.state()[b].q,solver.state()[b].qd,Vec3{a.phi,a.alpha,a.speed},Vec3{a.axial,a.tangential,a.coefficients.cl},Vec3{a.coefficients.cd,a.coefficients.cm,a.motion.position[0]},Vec3{a.motion.position[1],a.motion.position[2],a.load.force[0]},Vec3{a.load.force[1],a.load.force[2],a.load.moment[0]}})for(double x:v)out<<','<<x;
        out<<','<<a.load.moment[1]<<','<<a.load.moment[2]<<'\n';}
        if(solver.time()+c.dt/2>=c.duration)break;solver.step();}
    std::cout<<"Completed "<<solver.time()<<" s in "<<solver.step_number()<<" steps\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
