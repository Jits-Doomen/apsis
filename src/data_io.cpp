#include "apsis/data_io.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
namespace apsis::data {
std::vector<orbit::State> readStateCsv(const std::string& path){
 std::ifstream in(path);if(!in)throw std::runtime_error("cannot open state CSV: "+path);
 std::string line;if(!std::getline(in,line))return{};std::vector<orbit::State> out;std::size_t ln=1;
 while(std::getline(in,line)){++ln;if(line.empty())continue;std::replace(line.begin(),line.end(),',',' ');std::istringstream ss(line);double t,x,y,z,vx,vy,vz;
  if(!(ss>>t>>x>>y>>z>>vx>>vy>>vz))throw std::runtime_error("invalid state CSV row "+std::to_string(ln));
  out.push_back({{x,y,z},{vx,vy,vz},t,orbit::Frame::Inertial});
 }return out;
}
void writeStateCsv(const std::string&path,const std::vector<orbit::State>&s){
 std::ofstream out(path);if(!out)throw std::runtime_error("cannot write state CSV: "+path);
 out<<"epoch_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s\n"<<std::setprecision(17);
 for(const auto&q:s)out<<q.epoch_seconds<<','<<q.r.x<<','<<q.r.y<<','<<q.r.z<<','<<q.v.x<<','<<q.v.y<<','<<q.v.z<<'\n';
}
void writeXYCsv(const std::string&path,const std::vector<std::pair<double,double>>&d,const std::string&x,const std::string&y){
 std::ofstream out(path);if(!out)throw std::runtime_error("cannot write CSV: "+path);out<<x<<','<<y<<'\n'<<std::setprecision(17);for(auto [a,b]:d)out<<a<<','<<b<<'\n';
}
}
