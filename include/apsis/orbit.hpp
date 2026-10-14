#pragma once
#include <array>
#include <functional>
#include <string>
#include <vector>
namespace apsis::orbit {
struct Vec3 {
 double x=0,y=0,z=0;
 Vec3 operator+(Vec3 b)const; Vec3 operator-(Vec3 b)const; Vec3 operator*(double s)const; Vec3 operator/(double s)const;
 double dot(Vec3 b)const; Vec3 cross(Vec3 b)const; double norm()const; Vec3 unit()const;
};
enum class Frame { Inertial, BodyFixed };
struct State { Vec3 r; Vec3 v; double epoch_seconds=0; Frame frame=Frame::Inertial; };
struct Body { double mu; double radius=0; double j2=0; double rotation_rate=0; };
struct Forces {
 bool two_body=true, j2=false, drag=false, solar_radiation_pressure=false;
 double mass_kg=1, drag_coefficient=2.2, area_m2=0, density_kg_m3=0;
 double srp_pressure_n_m2=4.56e-6, reflectivity=1;
 Vec3 third_body_position_m{}; double third_body_mu=0;
 Vec3 thrust_newtons{}; // inertial force vector; divided by mass
};
struct PropagationOptions { double step_seconds=10; double max_step_seconds=60; };
Vec3 acceleration(const State& state,const Body& body,const Forces& forces);
State propagateRK4(const State& initial,const Body& body,const Forces& forces,
                   double duration_seconds,const PropagationOptions& options={});
State rotateZ(const State& state,double angle_radians,Frame target);
struct OrbitalElements { double semi_major_axis, eccentricity, inclination, raan, arg_periapsis, true_anomaly; };
OrbitalElements elementsFromState(const State& state,double mu);
State stateFromElements(const OrbitalElements& e,double mu,double epoch_seconds=0);
double specificEnergy(const State& s,double mu);
Vec3 angularMomentum(const State& s);
}
