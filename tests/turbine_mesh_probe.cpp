#include "turbine/mesh.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace turbine;
int main(int argc,char** argv){try{
    if(argc!=3)throw std::runtime_error("Usage: turbine_mesh_probe FIXTURE OUTPUT.csv");
    std::ofstream fixture(argv[1]),out(argv[2]);if(!fixture||!out)throw std::runtime_error("Cannot create test outputs");
    fixture<<std::setprecision(17);out<<std::setprecision(17)<<"kind,sample,node";for(int i=0;i<18;++i)out<<",v"<<i;out<<'\n';
    auto write=[&](std::ostream& s,const Vec3& v){for(double x:v)s<<x<<' ';};
    std::vector<ReferenceNode> sr,dr;std::vector<Vec3> sp,dp;
    for(int j=0;j<19;++j){double z=j==0?0:j==18?63:(j-.5)*63/17;
        sr.push_back({{0,0,z},euler_matrix({0,0,-.3*(1-z/63)})});if(j>0&&j<18)sp.push_back(sr.back().position);}
    for(int j=0;j<30;++j){double z=63.0*j/29;dr.push_back({{2*std::pow(z/63,2),0,z},euler_matrix({0,.02,-.25*(1-z/63)})});dp.push_back(dr.back().position);}
    fixture<<sr.size()<<' '<<dr.size()<<' '<<100<<'\n';
    for(const auto& nodes:{sr,dr})for(const auto& node:nodes){write(fixture,node.position);for(const auto& row:node.orientation)write(fixture,row);fixture<<'\n';}
    MotionMap motions(sr,dr);LoadMap loads(dp,sp);
    for(int n=0;n<100;++n){double t=n*.08;std::vector<Motion> sm(sr.size());std::vector<Vec3> sc,dc;std::vector<PointLoad> distributed;
        for(std::size_t j=0;j<sr.size();++j){const double z=sr[j].position[2];auto& m=sm[j];
            const auto r=euler_matrix({t,.03*std::sin(t+z/63),.01*z});m.orientation=multiply(sr[j].orientation,r);
            m.position=multiply(transpose(r),sr[j].position)+Vec3{.1*std::sin(t+z),.2*std::cos(t+z),.02*z};
            m.velocity={.1*z,std::sin(t+j),std::cos(t+j)};m.angular_velocity={1,.03*std::cos(t+z/63),.01};
            write(fixture,m.position);for(const auto& row:m.orientation)write(fixture,row);write(fixture,m.velocity);write(fixture,m.angular_velocity);fixture<<'\n';
            if(j>0&&j+1<sr.size())sc.push_back(m.position);
        }
        const auto dm=motions.transfer(sm);
        for(std::size_t j=0;j<dm.size();++j){dc.push_back(dm[j].position);PointLoad f{{1000+10.0*j,200*std::sin(t+j),20*std::cos(t+j)},{10*std::sin(t),5*std::cos(t),2}};
            distributed.push_back(f);write(fixture,f.force);write(fixture,f.moment);fixture<<'\n';
            out<<0<<','<<n<<','<<j;for(const auto& v:{dm[j].position,dm[j].velocity,dm[j].orientation[0],dm[j].orientation[1],dm[j].orientation[2],dm[j].angular_velocity})for(double x:v)out<<','<<x;out<<'\n';}
        const auto result=loads.transfer(distributed,dc,sc);
        for(std::size_t j=0;j<result.size();++j){out<<1<<','<<n<<','<<j;for(const auto& v:{result[j].force,result[j].moment})for(double x:v)out<<','<<x;for(int k=0;k<12;++k)out<<",0";out<<'\n';}
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
