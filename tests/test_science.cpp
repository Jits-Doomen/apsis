#include "apsis/numerics.hpp"
#include "apsis/linear_algebra.hpp"
#include "apsis/orbit.hpp"
#include "apsis/units.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

static void check(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
int main(){
 using namespace apsis;
 auto root=numerics::bisection([](double x){return x*x-2;},0,2);
 check(std::abs(root.root-std::sqrt(2.0))<1e-10,"bisection");
 check(std::abs(numerics::integrateAdaptiveSimpson([](double x){return x*x;},0,1)-1.0/3)<1e-10,"Simpson");
 auto ode=numerics::integrateRK45([](double,const numerics::State& y){return numerics::State{y[0]};},0,{1},1);
 check(std::abs(ode.y[0]-std::exp(1.0))<1e-7,"RK45");
 linalg::Matrix m(2,2);m(0,0)=3;m(0,1)=2;m(1,0)=1;m(1,1)=2;
 auto x=m.solve({5,5});check(std::abs(x[0])<1e-12&&std::abs(x[1]-2.5)<1e-12,"linear solve");
 check(std::abs(units::convert(1,"km","m")-1000)<1e-12,"unit conversion");
 orbit::Body earth{3.986004418e14,6378137,1.08262668e-3,7.292115e-5};
 orbit::State s{{7000000,0,0},{0,std::sqrt(earth.mu/7000000),0},0,orbit::Frame::Inertial};
 auto s2=orbit::propagateRK4(s,earth,{},60,{1,1});
 check(std::abs(s2.r.norm()-7000000)<1.0,"circular orbit radius drift");
 auto elems=orbit::elementsFromState(s,earth.mu);
 check(std::abs(elems.semi_major_axis-7000000)<1e-5,"orbital elements");
 std::cout<<"Apsis scientific tests passed\n";
}
