#pragma once

namespace ApsisPackages::Propagator {
    struct State3D {
        double x, y, z;
        double vx, vy, vz;
    };

    double compute_density(double alt);
    void gravitational_acceleration_spherical_harmonics(const State3D& s, double& ax, double& ay, double& az);
    void aerodynamic_drag_acceleration(const State3D& s, double cd_area_mass, double& ax, double& ay, double& az);
    void evaluate_derivatives(const State3D& s, double cd_a_m, State3D& ds);
    State3D propagate_rk4_step(const State3D& current, double dt, double cd_a_m);
}