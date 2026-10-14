#include "apsis/orbit.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace apsis::orbit {
namespace { constexpr double PI=3.14159265358979323846; }
Vec3 Vec3::operator+(Vec3 b)const{return{x+b.x,y+b.y,z+b.z};}
Vec3 Vec3::operator-(Vec3 b)const{return{x-b.x,y-b.y,z-b.z};}
Vec3 Vec3::operator*(double s)const{return{x*s,y*s,z*s};}
Vec3 Vec3::operator/(double s)const{if(s==0)throw std::domain_error("vector division by zero");return{x/s,y/s,z/s};}
double Vec3::dot(Vec3 b)const{return x*b.x+y*b.y+z*b.z;}
Vec3 Vec3::cross(Vec3 b)const{return{y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};}
double Vec3::norm()const{return std::sqrt(dot(*this));}
Vec3 Vec3::unit()const{double n=norm();if(!(n>0)||!std::isfinite(n))throw std::domain_error("cannot normalize zero or non-finite vector");return *this/n;}
Vec3 angularMomentum(const State&s){return s.r.cross(s.v);}
double specificEnergy(const State&s,double mu){double r=s.r.norm();if(r<=0||mu<=0)throw std::invalid_argument("invalid state or gravitational parameter");return .5*s.v.dot(s.v)-mu/r;}
Vec3 acceleration(const State&s,const Body&b,const Forces&f){
 const double r=s.r.norm(); if(!(r>0)||!(b.mu>0))throw std::invalid_argument("radius and mu must be positive");
 Vec3 a{};
 if(f.two_body)a=a+s.r*(-b.mu/(r*r*r));
 if(f.j2 && b.j2!=0 && b.radius>0){
   const double z2=s.r.z*s.r.z, r2=r*r, k=1.5*b.j2*b.mu*b.radius*b.radius/std::pow(r,5);
   a.x+=k*s.r.x*(5*z2/r2-1); a.y+=k*s.r.y*(5*z2/r2-1); a.z+=k*s.r.z*(5*z2/r2-3);
 }
 if(f.third_body_mu>0){
   Vec3 rb=f.third_body_position_m, d=rb-s.r; double dn=d.norm(), rn=rb.norm();
   if(dn<=0||rn<=0)throw std::invalid_argument("invalid third-body geometry");
   a=a+(d*(f.third_body_mu/(dn*dn*dn))-rb*(f.third_body_mu/(rn*rn*rn)));
 }
 if(f.drag){
   if(f.mass_kg<=0||f.area_m2<0||f.density_kg_m3<0)throw std::invalid_argument("invalid drag parameters");
   Vec3 vrel=s.v; if(b.rotation_rate!=0)vrel=vrel-Vec3{-b.rotation_rate*s.r.y,b.rotation_rate*s.r.x,0};
   const double speed=vrel.norm(); if(speed>0)a=a+vrel*(-.5*f.density_kg_m3*f.drag_coefficient*f.area_m2/f.mass_kg*speed);
 }
 if(f.solar_radiation_pressure){
   if(f.mass_kg<=0||f.area_m2<0)throw std::invalid_argument("invalid SRP parameters");
   const double scale=f.srp_pressure_n_m2*f.reflectivity*f.area_m2/f.mass_kg;
   // Simplified fixed-direction model: +X away from the Sun; caller must define frame accordingly.
   a.x+=scale;
 }
 if(f.mass_kg<=0) { if(f.thrust_newtons.norm()>0)throw std::invalid_argument("mass must be positive for thrust"); }
 else a=a+f.thrust_newtons/f.mass_kg;
 return a;
}
State propagateRK4(const State& initial,const Body& body,const Forces& forces,double duration,const PropagationOptions&o){
 if(duration<0||o.step_seconds<=0||o.max_step_seconds<=0)throw std::invalid_argument("invalid propagation duration/step");
 State s=initial;double elapsed=0;std::size_t n=0;
 while(elapsed<duration){
  double h=std::min({o.step_seconds,o.max_step_seconds,duration-elapsed});
  auto acc=[&](Vec3 r,Vec3 v){State q{r,v,s.epoch_seconds+elapsed,s.frame};return acceleration(q,body,forces);};
  Vec3 k1r=s.v,k1v=acc(s.r,s.v);
  Vec3 k2r=s.v+k1v*(h/2),k2v=acc(s.r+k1r*(h/2),k2r);
  Vec3 k3r=s.v+k2v*(h/2),k3v=acc(s.r+k2r*(h/2),k3r);
  Vec3 k4r=s.v+k3v*h,k4v=acc(s.r+k3r*h,k4r);
  s.r=s.r+(k1r+k2r*2+k3r*2+k4r)*(h/6);
  s.v=s.v+(k1v+k2v*2+k3v*2+k4v)*(h/6);
  elapsed+=h;s.epoch_seconds+=h;
  if(++n>10000000)throw std::runtime_error("propagation step budget exceeded");
 }
 return s;
}
State rotateZ(const State&s,double a,Frame target){
 if(s.frame==target)return s;double c=std::cos(a),sn=std::sin(a);
 auto rot=[&](Vec3 p){return Vec3{c*p.x-sn*p.y,sn*p.x+c*p.y,p.z};};
 State out=s;out.r=rot(s.r);out.v=rot(s.v);out.frame=target;return out;
}
OrbitalElements elementsFromState(const State&s,double mu){
 if(mu<=0)throw std::invalid_argument("mu must be positive");
 Vec3 r=s.r,v=s.v,h=r.cross(v);double rn=r.norm(),hn=h.norm();if(rn==0||hn==0)throw std::invalid_argument("degenerate state");
 Vec3 n=Vec3{0,0,1}.cross(h);Vec3 evec=(v.cross(h))/mu-r/rn;double ecc=evec.norm();
 double energy=.5*v.dot(v)-mu/rn;double a=std::abs(energy)<1e-15?INFINITY:-mu/(2*energy);
 double inc=std::acos(std::clamp(h.z/hn,-1.0,1.0));double raan=n.norm()>1e-12?std::atan2(n.y,n.x):0;
 double argp=0,nu=0;
 if(n.norm()>1e-12&&ecc>1e-12){argp=std::acos(std::clamp(n.dot(evec)/(n.norm()*ecc),-1.0,1.0));if(evec.z<0)argp=2*PI-argp;}
 if(ecc>1e-12){nu=std::acos(std::clamp(evec.dot(r)/(ecc*rn),-1.0,1.0));if(r.dot(v)<0)nu=2*PI-nu;}
 else if(n.norm()>1e-12){nu=std::acos(std::clamp(n.dot(r)/(n.norm()*rn),-1.0,1.0));if(r.z<0)nu=2*PI-nu;}
 return{a,ecc,inc,raan,argp,nu};
}
State stateFromElements(const OrbitalElements&e,double mu,double epoch){
 if(mu<=0||e.semi_major_axis<=0||e.eccentricity<0||e.eccentricity>=1)throw std::invalid_argument("stateFromElements currently requires an elliptic orbit");
 double p=e.semi_major_axis*(1-e.eccentricity*e.eccentricity),nu=e.true_anomaly;
 Vec3 rp{p*std::cos(nu)/(1+e.eccentricity*std::cos(nu)),p*std::sin(nu)/(1+e.eccentricity*std::cos(nu)),0};
 double q=std::sqrt(mu/p);Vec3 vp{-q*std::sin(nu),q*(e.eccentricity+std::cos(nu)),0};
 auto transform=[&](Vec3 u){double O=e.raan,w=e.arg_periapsis,i=e.inclination;
  double x=(std::cos(O)*std::cos(w)-std::sin(O)*std::sin(w)*std::cos(i))*u.x+(-std::cos(O)*std::sin(w)-std::sin(O)*std::cos(w)*std::cos(i))*u.y;
  double y=(std::sin(O)*std::cos(w)+std::cos(O)*std::sin(w)*std::cos(i))*u.x+(-std::sin(O)*std::sin(w)+std::cos(O)*std::cos(w)*std::cos(i))*u.y;
  double z=std::sin(w)*std::sin(i)*u.x+std::cos(w)*std::sin(i)*u.y;return Vec3{x,y,z};};
 return{transform(rp),transform(vp),epoch,Frame::Inertial};
}
}
