#pragma once

namespace ApsisPackages::Trajectory {
    double stumpff_c2(double z);
    double stumpff_c3(double z);
    void solve_lambert(
        double r1x, double r1y, double r1z,
        double r2x, double r2y, double r2z,
        double dt, int prograde,
        double& v1x, double& v1y, double& v1z,
        double& v2x, double& v2y, double& v2z
    );
}