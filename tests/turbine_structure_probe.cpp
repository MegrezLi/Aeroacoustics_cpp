#include "turbine/structure.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv){
    try {
        if(argc!=3)throw std::runtime_error("Usage: turbine_structure_probe CASE.fst OUTPUT.csv");
        turbine::Case c(argv[1]);turbine::BladeStructure structure(c);
        std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot create output");
        out<<"sample,blade,node,ax1,ax2,ax3,x,y,z,vx,vy,vz,nx,ny,nz\n"<<std::setprecision(17);
        for(int n=0;n<100;++n){double t=.037*n;for(int b=0;b<3;++b){const int k=b+1;
            turbine::ModalState s{{.8*std::sin(t+k),.1*std::cos(2*t+k),.4*std::sin(.7*t-k)},{.8*std::cos(t+k),-.2*std::sin(2*t+k),.28*std::cos(.7*t-k)}};
            std::vector<turbine::PointLoad> loads(structure.nodes.size());
            for(std::size_t j=1;j+1<loads.size();++j){const double w=structure.nodes[j].width;
                loads[j].force={w*(1000+10*j),w*200*std::sin(t+j),w*20*std::cos(t+j)};
                loads[j].moment={w*10*std::sin(t),w*5*std::cos(t),w*2};}
            const auto acceleration=structure.acceleration(t,b,s,loads);
            for(std::size_t j=0;j<structure.nodes.size();++j){const auto motion=structure.motion(t,b,s,structure.nodes[j]);out<<n<<','<<b<<','<<j;
                for(const auto& v:{acceleration,motion.position,motion.velocity,motion.orientation[2]})for(double x:v)out<<','<<x;out<<'\n';}
        }}
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
