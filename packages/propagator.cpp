#include "apsis/propagator.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>

namespace ApsisPackages::Propagator {

    constexpr double MU_EARTH = 3.986004418e14;
    constexpr double R_EARTH = 6378137.0;
    constexpr double J2 = 1.08262668e-3;
    constexpr double J3 = -2.53265e-6;
    constexpr double J4 = -1.61962e-6;
    constexpr double OMEGA_EARTH = 7.292115e-5;


    double compute_density(double alt) {
        if (alt < 0.0) alt = 0.0;
        if (alt > 1000000.0) return 0.0;
        double rho0 = 1.225;
        double H = 7200.0 + (alt / 100000.0) * 500.0;
        return rho0 * std::exp(-alt / H);
    }

    void gravitational_acceleration_spherical_harmonics(const State3D& s, double& ax, double& ay, double& az) {
        double r2 = s.x * s.x + s.y * s.y + s.z * s.z;
        double r = std::sqrt(r2);
        if (r == 0.0) throw std::domain_error("Singularity: center of Earth reached");

        double r_inv = 1.0 / r;
        double r2_inv = r_inv * r_inv;
        double r3_inv = r2_inv * r_inv * r;
        double z_r = s.z * r_inv;
        double z2_r2 = z_r * z_r;

        double f_two_body = -MU_EARTH * r3_inv;
        ax = f_two_body * s.x;
        ay = f_two_body * s.y;
        az = f_two_body * s.z;

        double r_eq_r = R_EARTH * r_inv;
        double r_eq_r2 = r_eq_r * r_eq_r;
        double j2_term = 1.5 * J2 * MU_EARTH * r3_inv * r_eq_r2;

        double five_z2 = 5.0 * z2_r2;
        ax += j2_term * s.x * (7.0 * z2_r2 - 1.0 - (five_z2 - 1.0));

        double term_j2 = 1.5 * J2 * MU_EARTH * std::pow(R_EARTH, 2.0) * std::pow(r, -5.0);
        ax += term_j2 * s.x * (5.0 * z2_r2 - 1.0);
        ay += term_j2 * s.y * (5.0 * z2_r2 - 1.0);
        az += term_j2 * s.z * (5.0 * z2_r2 - 3.0);

        double term_j3 = 2.5 * J3 * MU_EARTH * std::pow(R_EARTH, 3.0) * std::pow(r, -7.0);
        ax += term_j3 * s.x * (7.0 * z2_r2 - 3.0) * z_r;
        ay += term_j3 * s.y * (7.0 * z2_r2 - 3.0) * z_r;
        az += term_j3 * (7.0 * z2_r2 - 6.0) * z_r;
    }

    void aerodynamic_drag_acceleration(const State3D& s, double cd_area_mass, double& ax, double& ay, double& az) {
        double r = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
        double alt = r - R_EARTH;
        double rho = compute_density(alt);
        if (rho <= 0.0 || cd_area_mass <= 0.0) {
            ax = 0.0; ay = 0.0; az = 0.0;
            return;
        }

        double vx_rel = s.vx + OMEGA_EARTH * s.y;
        double vy_rel = s.vy - OMEGA_EARTH * s.x;
        double vz_rel = s.vz;
        double v_rel_mag = std::sqrt(vx_rel * vx_rel + vy_rel * vy_rel + vz_rel * vz_rel);

        double drag_coeff = -0.5 * rho * cd_area_mass * v_rel_mag;
        ax = drag_coeff * vx_rel;
        ay = drag_coeff * vy_rel;
        az = drag_coeff * vz_rel;
    }

    void evaluate_derivatives(const State3D& s, double cd_a_m, State3D& ds) {
        double gx, gy, gz;
        gravitational_acceleration_spherical_harmonics(s, gx, gy, gz);

        double dx_drag, dy_drag, dz_drag;
        aerodynamic_drag_acceleration(s, cd_a_m, dx_drag, dy_drag, dz_drag);

        ds.x = s.vx;
        ds.y = s.vy;
        ds.z = s.vz;
        ds.vx = gx + dx_drag;
        ds.vy = gy + dy_drag;
        ds.vz = gz + dz_drag;
    }

    State3D propagate_rk4_step(const State3D& current, double dt, double cd_a_m) {
        State3D k1, k2, k3, k4, temp;

        evaluate_derivatives(current, cd_a_m, k1);

        temp.x = current.x + 0.5 * dt * k1.x;
        temp.y = current.y + 0.5 * dt * k1.y;
        temp.z = current.z + 0.5 * dt * k1.z;
        temp.vx = current.vx + 0.5 * dt * k1.vx;
        temp.vy = current.vy + 0.5 * dt * k1.vy;
        temp.vz = current.vz + 0.5 * dt * k1.vz;
        evaluate_derivatives(temp, cd_a_m, k2);

        temp.x = current.x + 0.5 * dt * k2.x;
        temp.y = current.y + 0.5 * dt * k2.y;
        temp.z = current.z + 0.5 * dt * k2.z;
        temp.vx = current.vx + 0.5 * dt * k2.vx;
        temp.vy = current.vy + 0.5 * dt * k2.vy;
        temp.vz = current.vz + 0.5 * dt * k2.vz;
        evaluate_derivatives(temp, cd_a_m, k3);

        temp.x = current.x + dt * k3.x;
        temp.y = current.y + dt * k3.y;
        temp.z = current.z + dt * k3.z;
        temp.vx = current.vx + dt * k3.vx;
        temp.vy = current.vy + dt * k3.vy;
        temp.vz = current.vz + dt * k3.vz;
        evaluate_derivatives(temp, cd_a_m, k4);

        State3D next_state;
        next_state.x = current.x + (dt / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x);
        next_state.y = current.y + (dt / 6.0) * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y);
        next_state.z = current.z + (dt / 6.0) * (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z);
        next_state.vx = current.vx + (dt / 6.0) * (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx);
        next_state.vy = current.vy + (dt / 6.0) * (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy);
        next_state.vz = current.vz + (dt / 6.0) * (k1.vz + 2.0 * k2.vz + 2.0 * k3.vz + k4.vz);

        return next_state;
    }

}