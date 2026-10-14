#include "apsis/ephemerides.hpp"
#include <cmath>
#include <stdexcept>
#include <array>

namespace ApsisPackages::Ephemerides {

    constexpr double PI = 3.1415926535897932384626433832795;
    constexpr double TWO_PI = 2.0 * PI;
    constexpr double DEG_TO_RAD = PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / PI;
    constexpr double AU_TO_M = 149597870700.0;

    double to_julian_date(int year, int month, int day, double hour, double minute, double second) {
        if (month <= 2) {
            year -= 1;
            month += 12;
        }
        int A = year / 100;
        int B = 2 - A + (A / 4);
        double day_fraction = (hour + (minute / 60.0) + (second / 3600.0)) / 24.0;
        return static_cast<double>(int(365.25 * (year + 4716))) +
            static_cast<double>(int(30.6001 * (month + 1))) +
            day_fraction + day + B - 1524.5;
    }

    double centuries_from_j2000(double jd) {
        return (jd - 2451545.0) / 36525.0;
    }

    double greenwich_mean_sidereal_time(double jd) {
        double T = centuries_from_j2000(jd);
        double gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) +
            0.000387933 * T * T - (T * T * T) / 38710000.0;
        gmst_deg = std::fmod(gmst_deg, 360.0);
        if (gmst_deg < 0.0) gmst_deg += 360.0;
        return gmst_deg * DEG_TO_RAD;
    }

    void get_sun_position_ecliptic(double jd, double& x_ecl, double& y_ecl, double& z_ecl) {
        double T = centuries_from_j2000(jd);
        double M = std::fmod(357.52911 + 35999.05029 * T, 360.0) * DEG_TO_RAD;
        double L0 = std::fmod(280.46646 + 36000.76983 * T, 360.0) * DEG_TO_RAD;

        double C_sun = (1.914602 - 0.004817 * T) * std::sin(M) +
            (0.019993 - 0.000101 * T) * std::sin(2.0 * M) +
            0.000289 * std::sin(3.0 * M);
        double lambda = L0 + C_sun * DEG_TO_RAD;

        double r_au = 1.00014010 - 0.01673116 * std::cos(M) - 0.0001395 * std::cos(2.0 * M);
        double r = r_au * AU_TO_M;

        x_ecl = r * std::cos(lambda);
        y_ecl = r * std::sin(lambda);
        z_ecl = 0.0;
    }

    void ecliptic_to_equatorial(double jd, double x_ecl, double y_ecl, double z_ecl, double& x_eq, double& y_eq, double& z_eq) {
        double T = centuries_from_j2000(jd);
        double eps = (23.43929111 - (46.8150 * T + 0.00059 * T * T - 0.001813 * T * T * T) / 3600.0) * DEG_TO_RAD;

        double coseps = std::cos(eps);
        double sineps = std::sin(eps);

        x_eq = x_ecl;
        y_eq = y_ecl * coseps - z_ecl * sineps;
        z_eq = y_ecl * sineps + z_ecl * coseps;
    }

    void eci_to_ecef(double jd, double x_eci, double y_eci, double z_eci, double& x_ecef, double& y_ecef, double& z_ecef) {
        double gmst = greenwich_mean_sidereal_time(jd);
        double cos_g = std::cos(gmst);
        double sin_g = std::sin(gmst);

        x_ecef = cos_g * x_eci + sin_g * y_eci;
        y_ecef = -sin_g * x_eci + cos_g * y_eci;
        z_ecef = z_eci;
    }

    void ecef_to_eci(double jd, double x_ecef, double y_ecef, double z_ecef, double& x_eci, double& y_eci, double& z_eci) {
        double gmst = greenwich_mean_sidereal_time(jd);
        double cos_g = std::cos(gmst);
        double sin_g = std::sin(gmst);

        x_eci = cos_g * x_ecef - sin_g * y_ecef;
        y_eci = sin_g * x_ecef + cos_g * y_ecef;
        z_eci = z_ecef;
    }

}