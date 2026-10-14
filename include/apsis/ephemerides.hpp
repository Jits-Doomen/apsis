#pragma once

namespace ApsisPackages::Ephemerides {
    double to_julian_date(int year, int month, int day, double hour, double minute, double second);
    double centuries_from_j2000(double jd);
    double greenwich_mean_sidereal_time(double jd);
    void get_sun_position_ecliptic(double jd, double& x_ecl, double& y_ecl, double& z_ecl);
    void ecliptic_to_equatorial(double jd, double x_ecl, double y_ecl, double z_ecl, double& x_eq, double& y_eq, double& z_eq);
    void eci_to_ecef(double jd, double x_eci, double y_eci, double z_eci, double& x_ecef, double& y_ecef, double& z_ecef);
    void ecef_to_eci(double jd, double x_ecef, double y_ecef, double z_ecef, double& x_eci, double& y_eci, double& z_eci);
}