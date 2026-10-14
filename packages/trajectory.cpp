#include "apsis/trajectory.hpp"
#include <cmath>
#include <stdexcept>

namespace ApsisPackages::Trajectory {

    constexpr double MU_EARTH = 3.986004418e14;
    constexpr double EPSILON = 1e-8;
    constexpr int MAX_ITERATIONS = 100;

    double stumpff_c2(double z) {
        if (z > EPSILON) {
            double sz = std::sqrt(z);
            return (1.0 - std::cos(sz)) / z;
        }
        else if (z < -EPSILON) {
            double sz = std::sqrt(-z);
            return (1.0 - std::cosh(sz)) / z;
        }
        else {
            return 0.5 - (z / 24.0) + (z * z / 720.0);
        }
    }

    double stumpff_c3(double z) {
        if (z > EPSILON) {
            double sz = std::sqrt(z);
            return (sz - std::sin(sz)) / (z * sz);
        }
        else if (z < -EPSILON) {
            double sz = std::sqrt(-z);
            return (std::sinh(sz) - sz) / (-z * sz);
        }
        else {
            return (1.0 / 6.0) - (z / 120.0) + (z * z / 5040.0);
        }
    }

    void solve_lambert(
        double r1x, double r1y, double r1z,
        double r2x, double r2y, double r2z,
        double dt, int prograde,
        double& v1x, double& v1y, double& v1z,
        double& v2x, double& v2y, double& v2z
    ) {
        double r1 = std::sqrt(r1x * r1x + r1y * r1y + r1z * r1z);
        double r2 = std::sqrt(r2x * r2x + r2y * r2y + r2z * r2z);

        double cross_z = r1x * r2y - r1y * r2x;
        double dot_val = (r1x * r2x + r1y * r2y + r1z * r2z) / (r1 * r2);
        if (dot_val > 1.0) dot_val = 1.0;
        if (dot_val < -1.0) dot_val = -1.0;

        double dtheta = std::acos(dot_val);
        if (prograde && cross_z < 0.0) dtheta = 2.0 * 3.141592653589793 - dtheta;
        if (!prograde && cross_z >= 0.0) dtheta = 2.0 * 3.141592653589793 - dtheta;

        double A = std::sin(dtheta) * std::sqrt((r1 * r2) / (1.0 - std::cos(dtheta)));
        if (std::abs(A) < EPSILON) {
            throw std::invalid_argument("Lambert error: transfer angle is 0 or 180 degrees");
        }

        double z = 0.0;
        double z_low = -100.0;
        double z_high = 100.0;
        double dt_calc = 0.0;

        for (int i = 0; i < MAX_ITERATIONS; ++i) {
            double c2 = stumpff_c2(z);
            double c3 = stumpff_c3(z);
            double y = r1 + r2 + (A * (z * c3 - 1.0)) / std::sqrt(c2);

            if (A > 0.0 && y < 0.0) {
                z_low = z;
                z += 10.0;
                continue;
            }

            dt_calc = (std::pow(y, 1.5) * c3 + A * std::sqrt(y)) / std::sqrt(MU_EARTH);

            if (std::abs(dt_calc - dt) < 1e-6) {
                break;
            }

            if (dt_calc <= dt) {
                z_low = z;
            }
            else {
                z_high = z;
            }
            z = 0.5 * (z_low + z_high);
        }

        double c2 = stumpff_c2(z);
        double c3 = stumpff_c3(z);
        double y = r1 + r2 + (A * (z * c3 - 1.0)) / std::sqrt(c2);

        double f = 1.0 - y / r1;
        double g = A * std::sqrt(y / MU_EARTH);
        double g_dot = 1.0 - y / r2;

        v1x = (r2x - f * r1x) / g;
        v1y = (r2y - f * r1y) / g;
        v1z = (r2z - f * r1z) / g;

        v2x = (g_dot * r2x - r1x) / g;
        v2y = (g_dot * r2y - r1y) / g;
        v2z = (g_dot * r2z - r1z) / g;
    }

}