#include "apsis/formulas.hpp"
#include <cmath>
#include <stdexcept>

namespace ApsisPackages::Formulas {
    constexpr double PI = 3.14159265358979323846264338327950288419716939937510582097494459230;

    // Speed of light in vacuum (m/s)
    constexpr double C = 2.99792458e8;

    // Gravitational constant (m^3 / (kg * s^2))
    constexpr double G = 6.67430e-11;

    // Planck constant (J * s)
    constexpr double H = 6.62607015e-34;

    // Reduced Planck constant (J * s)
    constexpr double H_BAR = 1.054571817e-34;

    // Stefan-Boltzmann constant (W / (m^2 * K^4))
    constexpr double SIGMA = 5.670374e-8;

    // Boltzmann constant (J / K)
    constexpr double K_B = 1.380649e-23;

    // Universal gas constant (J / (mol * K))
    constexpr double R_GAS = 8.31446261815324;

    // Standard gravity (m/s^2)
    constexpr double G0 = 9.80665;

    // Permittivity of free space (F / m)
    constexpr double EPSILON_0 = 8.8541878128e-12;

    // Permeability of free space (N / A^2)
    constexpr double MU_0 = 1.25663706e-6;

    // Coulomb constant (N * m^2 / C^2)
    constexpr double K_E = 8.9875517923e9;

    // Wien's displacement constant (m * K)
    constexpr double WIEN_B = 2.897771955e-3;

    // Proton mass (kg)
    constexpr double M_P = 1.67262192369e-27;

    // Hydrogen atom mass (kg)
    constexpr double M_H = 1.6735575e-27;

    // Thomson cross-section (m^2)
    constexpr double SIGMA_T = 6.6524587321e-29;


    double projectile_max_heigth(double v0, double theta, double g) {
        if (g <= 0.0) {
            throw std::invalid_argument("projectile_max_heigth: g must be positive");
        }
        const double sin_theta = std::sin(theta);
        return (v0 * v0 * sin_theta * sin_theta) / (2.0 * g);
    }

    double projectile_range(double v0, double theta, double g) {
        if (g <= 0.0) {
            throw std::invalid_argument("projectile_range: g must be positive");
        }
        return (v0 * v0 * std::sin(2.0 * theta)) / g;
    }

    double projectile_time_of_flight(double v0, double theta, double g) {
        if (g <= 0.0) {
            throw std::invalid_argument("projectile_time_of_flight: g must be positive");
        }
        return (2.0 * v0 * std::sin(theta)) / g;
    }

    double average_velocity(double d, double t) {
        if (t == 0.0) {
            throw std::invalid_argument("average_velocity: time must not be zero");
        }
        return d / t;
    }

    double average_acceleration(double delta_v, double delta_t) {
        if (delta_t == 0.0) {
            throw std::invalid_argument("average_acceleration: time interval must not be zero");
        }
        return delta_v / delta_t;
    }

    double kinematic_final_velocity(double v0, double a, double t) {
        return v0 + a * t;
    }

    double kinematic_position(double x0, double v0, double a, double t) {
        return x0 + v0 * t + 0.5 * a * t * t;
    }

    double kinematic_velocity_squared(double v0, double a, double delta_x) {
        return v0 * v0 + 2.0 * a * delta_x;
    }

    double newton_second_law(double m, double a) {
        return m * a;
    }

    double weight(double m, double g) {
        return m * g;
    }

    double gravitational_force(double m1, double m2, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("gravitational_force: distance must be positive");
        }
        return (G * m1 * m2) / (r * r);
    }

    double linear_momentum(double m, double v) {
        return m * v;
    }

    double impulse(double F, double delta_t) {
        return F * delta_t;
    }

    double friction_force(double mu, double N) {
        return mu * N;
    }

    double work_done(double F, double d, double theta) {
        return F * d * std::cos(theta);
    }

    double kinetic_energy(double m, double v) {
        return 0.5 * m * v * v;
    }

    double gravitational_potential_energy(double m, double g, double h) {
        return m * g * h;
    }

    double mechanical_power(double W, double t) {
        if (t <= 0.0) {
            throw std::invalid_argument("mechanical_power: time must be positive");
        }
        return W / t;
    }

    double hookes_law(double k, double x) {
        return -k * x;
    }

    double elastic_potential_energy(double k, double x) {
        return 0.5 * k * x * x;
    }

    double angular_velocity(double delta_theta, double delta_t) {
        if (delta_t == 0.0) {
            throw std::invalid_argument("angular_velocity: time interval must not be zero");
        }
        return delta_theta / delta_t;
    }

    double angular_acceleration(double delta_omega, double delta_t) {
        if (delta_t == 0.0) {
            throw std::invalid_argument("angular_acceleration: time interval must not be zero");
        }
        return delta_omega / delta_t;
    }

    double tangential_velocity(double omega, double r) {
        return omega * r;
    }

    double centripetal_acceleration(double v, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("centripetal_acceleration: radius must be positive");
        }
        return (v * v) / r;
    }

    double centripetal_force(double m, double v, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("centripetal_force: radius must be positive");
        }
        return (m * v * v) / r;
    }

    double torque(double r, double F, double theta) {
        return r * F * std::sin(theta);
    }

    double moment_of_inertia_point_mass(double m, double r) {
        return m * r * r;
    }

    double moment_of_inertia_solid_disk(double m, double r) {
        return 0.5 * m * r * r;
    }

    double moment_of_inertia_solid_sphere(double m, double r) {
        return (2.0 / 5.0) * m * r * r;
    }

    double moment_of_inertia_thin_rod_center(double m, double L) {
        return (1.0 / 12.0) * m * L * L;
    }

    double angular_momentum(double I, double omega) {
        return I * omega;
    }

    double rotational_kinetic_energy(double I, double omega) {
        return 0.5 * I * omega * omega;
    }

    double parallel_axis_theorem(double I_cm, double m, double d) {
        return I_cm + m * d * d;
    }

    double shm_spring_period(double m, double k) {
        if (m <= 0.0 || k <= 0.0) {
            throw std::invalid_argument("shm_spring_period: mass and spring constant must be positive");
        }
        return 2.0 * PI * std::sqrt(m / k);
    }

    double simple_pendulum_period(double L, double g) {
        if (L <= 0.0 || g <= 0.0) {
            throw std::invalid_argument("simple_pendulum_period: length and gravity must be positive");
        }
        return 2.0 * PI * std::sqrt(L / g);
    }

    double shm_position(double A, double omega, double t, double phi) {
        return A * std::cos(omega * t + phi);
    }

    double shm_angular_frequency(double k, double m) {
        if (k <= 0.0 || m <= 0.0) {
            throw std::invalid_argument("shm_angular_frequency: spring constant and mass must be positive");
        }
        return std::sqrt(k / m);
    }

    double damped_oscillation_amplitude(double A0, double b, double m, double t) {
        if (m <= 0.0) {
            throw std::invalid_argument("damped_oscillation_amplitude: mass must be positive");
        }
        return A0 * std::exp((-b * t) / (2.0 * m));
    }

    double ideal_gas_pressure(double n, double T, double V) {
        if (V <= 0.0) {
            throw std::invalid_argument("ideal_gas_pressure: volume must be positive");
        }
        return (n * R_GAS * T) / V;
    }

    double first_law_thermodynamics(double Q, double W) {
        return Q - W;
    }

    double sensible_heat(double m, double c, double delta_T) {
        return m * c * delta_T;
    }

    double latent_heat(double m, double L) {
        return m * L;
    }

    double carnot_efficiency(double T_c, double T_h) {
        if (T_h <= 0.0 || T_c < 0.0 || T_c >= T_h) {
            throw std::invalid_argument("carnot_efficiency: invalid temperature parameters");
        }
        return 1.0 - (T_c / T_h);
    }

    double entropy_change_reversible(double Q, double T) {
        if (T <= 0.0) {
            throw std::invalid_argument("entropy_change_reversible: absolute temperature must be positive");
        }
        return Q / T;
    }

    double statistical_entropy(double W) {
        if (W <= 0.0) {
            throw std::invalid_argument("statistical_entropy: microstates W must be positive");
        }
        return K_B * std::log(W);
    }

    double stefan_boltzmann_power(double A, double T) {
        if (A < 0.0 || T < 0.0) {
            throw std::invalid_argument("stefan_boltzmann_power: area and temperature cannot be negative");
        }
        return SIGMA * A * std::pow(T, 4.0);
    }

    double wiens_displacement_law(double T) {
        if (T <= 0.0) {
            throw std::invalid_argument("wiens_displacement_law: temperature must be positive");
        }
        return WIEN_B / T;
    }

    double linear_thermal_expansion(double alpha, double L0, double delta_T) {
        return alpha * L0 * delta_T;
    }

    double fourier_conduction_rate(double k, double A, double delta_T, double d) {
        if (d <= 0.0) {
            throw std::invalid_argument("fourier_conduction_rate: thickness must be positive");
        }
        return (k * A * delta_T) / d;
    }

    double rms_molecular_speed(double T, double M_molar) {
        if (T < 0.0 || M_molar <= 0.0) {
            throw std::invalid_argument("rms_molecular_speed: temperature cannot be negative and molar mass must be positive");
        }
        return std::sqrt((3.0 * R_GAS * T) / M_molar);
    }

    double isobaric_work(double P, double delta_V) {
        return P * delta_V;
    }

    double coulomb_force(double q1, double q2, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("coulomb_force: separation distance must be positive");
        }
        return (K_E * q1 * q2) / (r * r);
    }

    double electric_field_point_charge(double Q_charge, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("electric_field_point_charge: distance must be positive");
        }
        return (K_E * Q_charge) / (r * r);
    }

    double electric_force_on_charge(double q, double E_field) {
        return q * E_field;
    }

    double electric_potential_point_charge(double Q_charge, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("electric_potential_point_charge: distance must be positive");
        }
        return (K_E * Q_charge) / r;
    }

    double parallel_plate_capacitance(double A, double d) {
        if (d <= 0.0) {
            throw std::invalid_argument("parallel_plate_capacitance: plate separation must be positive");
        }
        return (EPSILON_0 * A) / d;
    }

    double capacitor_stored_energy(double C_cap, double V) {
        return 0.5 * C_cap * V * V;
    }

    double ohms_law_voltage(double I, double R) {
        return I * R;
    }

    double electrical_power(double V, double I) {
        return V * I;
    }

    double conductor_resistance(double rho, double L, double A) {
        if (A <= 0.0) {
            throw std::invalid_argument("conductor_resistance: area must be positive");
        }
        return (rho * L) / A;
    }

    double amperes_law_magnetic_field_circulation(double I_enc) {
        return MU_0 * I_enc;
    }

    double faradays_law_emf(double delta_phi_B, double delta_t) {
        if (delta_t <= 0.0) {
            throw std::invalid_argument("faradays_law_emf: time interval must be positive");
        }
        return -delta_phi_B / delta_t;
    }

    double self_inductance_emf(double L_ind, double delta_I, double delta_t) {
        if (delta_t <= 0.0) {
            throw std::invalid_argument("self_inductance_emf: time interval must be positive");
        }
        return -L_ind * (delta_I / delta_t);
    }

    double inductor_stored_energy(double L_ind, double I) {
        return 0.5 * L_ind * I * I;
    }

    double speed_of_em_waves() {
        return 1.0 / std::sqrt(MU_0 * EPSILON_0);
    }

    double wave_speed(double f, double lambda) {
        return f * lambda;
    }

    double wave_period(double f) {
        if (f <= 0.0) {
            throw std::invalid_argument("wave_period: frequency must be positive");
        }
        return 1.0 / f;
    }

    double snells_law_angle2(double n1, double theta1, double n2) {
        if (n2 <= 0.0) {
            throw std::invalid_argument("snells_law_angle2: refractive index n2 must be positive");
        }
        const double sin_theta2 = (n1 * std::sin(theta1)) / n2;
        if (sin_theta2 < -1.0 || sin_theta2 > 1.0) {
            throw std::domain_error("snells_law_angle2: total internal reflection occurs");
        }
        return std::asin(sin_theta2);
    }

    double critical_angle(double n1, double n2) {
        if (n1 <= 0.0 || n2 <= 0.0 || n2 > n1) {
            throw std::invalid_argument("critical_angle: invalid refractive indices (requires n1 >= n2 > 0)");
        }
        return std::asin(n2 / n1);
    }

    double doppler_effect_observed_frequency(double f_source, double v_wave, double v_observer, double v_source) {
        const double denom = v_wave - v_source;
        if (denom == 0.0) {
            throw std::invalid_argument("doppler_effect_observed_frequency: denominator evaluated to zero");
        }
        return f_source * ((v_wave + v_observer) / denom);
    }

    double thin_lens_image_distance(double f, double d_o) {
        if (d_o == f) {
            throw std::invalid_argument("thin_lens_image_distance: object distance equal to focal length");
        }
        return (f * d_o) / (d_o - f);
    }

    double lens_magnification(double d_i, double d_o) {
        if (d_o == 0.0) {
            throw std::invalid_argument("lens_magnification: object distance cannot be zero");
        }
        return -d_i / d_o;
    }

    double wave_intensity(double P, double A) {
        if (A <= 0.0) {
            throw std::invalid_argument("wave_intensity: area must be positive");
        }
        return P / A;
    }

    double sound_level_decibels(double I) {
        constexpr double I0 = 1.0e-12;
        if (I <= 0.0) {
            throw std::invalid_argument("sound_level_decibels: intensity must be positive");
        }
        return 10.0 * std::log10(I / I0);
    }

    double mass_energy_equivalence(double m) {
        return m * C * C;
    }

    double lorentz_factor(double v) {
        if (std::abs(v) >= C) {
            throw std::invalid_argument("lorentz_factor: velocity must be less than c");
        }
        return 1.0 / std::sqrt(1.0 - (v * v) / (C * C));
    }

    double time_dilation(double delta_t0, double gamma) {
        return gamma * delta_t0;
    }

    double length_contraction(double L0, double gamma) {
        if (gamma <= 0.0) {
            throw std::invalid_argument("length_contraction: gamma must be positive");
        }
        return L0 / gamma;
    }

    double relativistic_momentum(double m, double v, double gamma) {
        return gamma * m * v;
    }

    double relativistic_total_energy(double m, double gamma) {
        return gamma * m * C * C;
    }

    double energy_momentum_relation(double p, double m) {
        const double pc = p * C;
        const double mc2 = m * C * C;
        return std::sqrt(pc * pc + mc2 * mc2);
    }

    double photon_energy_frequency(double f) {
        return H * f;
    }

    double de_broglie_wavelength(double p) {
        if (p <= 0.0) {
            throw std::invalid_argument("de_broglie_wavelength: momentum must be positive");
        }
        return H / p;
    }

    double radioactive_decay_remaining(double N0, double lambda_decay, double t) {
        return N0 * std::exp(-lambda_decay * t);
    }

    double radioactive_half_life(double lambda_decay) {
        if (lambda_decay <= 0.0) {
            throw std::invalid_argument("radioactive_half_life: decay constant must be positive");
        }
        return std::log(2.0) / lambda_decay;
    }

    double keplers_third_law(double a, double M) {
        if (a <= 0.0) {
            throw std::invalid_argument("keplers_third_law: semi-major axis must be positive");
        }
        if (M <= 0.0) {
            throw std::invalid_argument("keplers_third_law: central mass must be positive");
        }
        const double numerator = 4.0 * PI * PI * a * a * a;
        const double denominator = G * M;
        return std::sqrt(numerator / denominator);
    }

    double stellar_luminosity(double R, double T) {
        if (R < 0.0) {
            throw std::invalid_argument("stellar_luminosity: radius cannot be negative");
        }
        if (T < 0.0) {
            throw std::invalid_argument("stellar_luminosity: temperature cannot be negative");
        }
        return 4.0 * PI * R * R * SIGMA * T * T * T * T;
    }

    double flux_luminosity_relation(double L, double d) {
        if (d <= 0.0) {
            throw std::invalid_argument("flux_luminosity_relation: distance must be positive");
        }
        return L / (4.0 * PI * d * d);
    }

    double distance_modulus(double d) {
        if (d <= 0.0) {
            throw std::invalid_argument("distance_modulus: distance must be positive");
        }
        return 5.0 * std::log10(d / 10.0);
    }

    double parallax_distance(double p_arcsec) {
        if (p_arcsec <= 0.0) {
            throw std::invalid_argument("parallax_distance: parallax angle must be positive");
        }
        return 1.0 / p_arcsec;
    }

    double hubbles_law_velocity(double H0, double d) {
        return H0 * d;
    }

    double cosmological_redshift(double lambda_obs, double lambda_rest) {
        if (lambda_rest <= 0.0) {
            throw std::invalid_argument("cosmological_redshift: rest wavelength must be positive");
        }
        return (lambda_obs - lambda_rest) / lambda_rest;
    }

    double schwarzschild_radius(double M) {
        if (M <= 0.0) {
            throw std::invalid_argument("schwarzschild_radius: mass must be positive");
        }
        return (2.0 * G * M) / (C * C);
    }

    double eddington_luminosity(double M) {
        if (M <= 0.0) {
            throw std::invalid_argument("eddington_luminosity: mass must be positive");
        }
        return (4.0 * PI * G * M * M_P * C) / SIGMA_T;
    }

    double vis_viva_velocity(double G_const, double M, double r, double a) {
        if (r <= 0.0 || a <= 0.0) {
            throw std::invalid_argument("vis_viva_velocity: radial position and semi-major axis must be positive");
        }
        return std::sqrt(G_const * M * ((2.0 / r) - (1.0 / a)));
    }

    double circular_orbital_velocity(double M, double r) {
        if (r <= 0.0 || M <= 0.0) {
            throw std::invalid_argument("circular_orbital_velocity: radius and mass must be positive");
        }
        return std::sqrt((G * M) / r);
    }

    double escape_velocity(double M, double r) {
        if (r <= 0.0 || M <= 0.0) {
            throw std::invalid_argument("escape_velocity: distance and mass must be positive");
        }
        return std::sqrt((2.0 * G * M) / r);
    }

    double specific_orbital_energy(double v, double M, double r) {
        if (r <= 0.0) {
            throw std::invalid_argument("specific_orbital_energy: distance must be positive");
        }
        return (v * v / 2.0) - ((G * M) / r);
    }

    double orbital_eccentricity_from_apsides(double r_apo, double r_peri) {
        if (r_apo + r_peri == 0.0) {
            throw std::invalid_argument("orbital_eccentricity_from_apsides: denominator sum cannot be zero");
        }
        return (r_apo - r_peri) / (r_apo + r_peri);
    }

    double semi_major_axis_from_apsides(double r_apo, double r_peri) {
        return (r_apo + r_peri) / 2.0;
    }

    double hohmann_transfer_delta_v1(double M, double r1, double r2) {
        if (r1 <= 0.0 || r2 <= 0.0) {
            throw std::invalid_argument("hohmann_transfer_delta_v1: orbital radii must be positive");
        }
        const double v_circ1 = std::sqrt((G * M) / r1);
        return v_circ1 * (std::sqrt((2.0 * r2) / (r1 + r2)) - 1.0);
    }

    double hohmann_transfer_delta_v2(double M, double r1, double r2) {
        if (r1 <= 0.0 || r2 <= 0.0) {
            throw std::invalid_argument("hohmann_transfer_delta_v2: orbital radii must be positive");
        }
        const double v_circ2 = std::sqrt((G * M) / r2);
        return v_circ2 * (1.0 - std::sqrt((2.0 * r1) / (r1 + r2)));
    }

    double sphere_of_influence_radius(double a, double m_small, double M_large) {
        if (a <= 0.0 || M_large <= 0.0) {
            throw std::invalid_argument("sphere_of_influence_radius: invalid orbital or mass inputs");
        }
        return a * std::pow(m_small / M_large, 0.4);
    }

    double tsiolkovsky_rocket_equation(double v_e, double m0, double mf) {
        if (mf <= 0.0 || m0 < mf) {
            throw std::invalid_argument("tsiolkovsky_rocket_equation: initial mass must be greater than dry mass");
        }
        return v_e * std::log(m0 / mf);
    }

    double effective_exhaust_velocity(double I_sp, double g0_acc) {
        return I_sp * g0_acc;
    }

    double specific_impulse(double F_thrust, double m_dot, double g0_acc) {
        if (m_dot <= 0.0 || g0_acc <= 0.0) {
            throw std::invalid_argument("specific_impulse: mass flow rate and g0 must be positive");
        }
        return F_thrust / (m_dot * g0_acc);
    }

    double thrust_equation(double m_dot, double v_e, double p_e, double p_a, double A_e) {
        return (m_dot * v_e) + ((p_e - p_a) * A_e);
    }

    double rocket_mass_ratio(double m0, double mf) {
        if (mf <= 0.0) {
            throw std::invalid_argument("rocket_mass_ratio: final mass must be positive");
        }
        return m0 / mf;
    }

    double propellant_mass_fraction(double m_p, double m0) {
        if (m0 <= 0.0) {
            throw std::invalid_argument("propellant_mass_fraction: initial mass must be positive");
        }
        return m_p / m0;
    }

    double thrust_to_weight_ratio(double F_thrust, double m, double g0_acc) {
        if (m <= 0.0 || g0_acc <= 0.0) {
            throw std::invalid_argument("thrust_to_weight_ratio: mass and g0 must be positive");
        }
        return F_thrust / (m * g0_acc);
    }

    double rocket_burn_time(double m_p, double m_dot) {
        if (m_dot <= 0.0) {
            throw std::invalid_argument("rocket_burn_time: mass flow rate must be positive");
        }
        return m_p / m_dot;
    }

    double normal_stress(double F, double A) {
        if (A <= 0.0) {
            throw std::invalid_argument("normal_stress: cross-sectional area must be positive");
        }
        return F / A;
    }

    double normal_strain(double delta_L, double L0) {
        if (L0 <= 0.0) {
            throw std::invalid_argument("normal_strain: original length must be positive");
        }
        return delta_L / L0;
    }

    double youngs_modulus(double sigma, double epsilon) {
        if (epsilon == 0.0) {
            throw std::invalid_argument("youngs_modulus: strain cannot be zero");
        }
        return sigma / epsilon;
    }

    double shear_modulus(double tau, double gamma_strain) {
        if (gamma_strain == 0.0) {
            throw std::invalid_argument("shear_modulus: shear strain cannot be zero");
        }
        return tau / gamma_strain;
    }

    double poissons_ratio(double epsilon_transverse, double epsilon_axial) {
        if (epsilon_axial == 0.0) {
            throw std::invalid_argument("poissons_ratio: axial strain cannot be zero");
        }
        return -epsilon_transverse / epsilon_axial;
    }

    double bending_stress(double M_moment, double y, double I) {
        if (I <= 0.0) {
            throw std::invalid_argument("bending_stress: second moment of area must be positive");
        }
        return (M_moment * y) / I;
    }

    double torsional_shear_stress(double T_torque, double r, double J) {
        if (J <= 0.0) {
            throw std::invalid_argument("torsional_shear_stress: polar moment of inertia must be positive");
        }
        return (T_torque * r) / J;
    }

    double angle_of_twist(double T_torque, double L, double J, double G_modulus) {
        if (J <= 0.0 || G_modulus <= 0.0) {
            throw std::invalid_argument("angle_of_twist: polar moment and shear modulus must be positive");
        }
        return (T_torque * L) / (J * G_modulus);
    }

    double euler_buckling_load(double E_modulus, double I, double K_factor, double L) {
        const double effective_len = K_factor * L;
        if (effective_len <= 0.0) {
            throw std::invalid_argument("euler_buckling_load: effective length must be positive");
        }
        return (PI * PI * E_modulus * I) / (effective_len * effective_len);
    }

    double factor_of_safety(double sigma_yield, double sigma_allowable) {
        if (sigma_allowable <= 0.0) {
            throw std::invalid_argument("factor_of_safety: allowable stress must be positive");
        }
        return sigma_yield / sigma_allowable;
    }

    double thermal_stress_constrained(double E_modulus, double alpha, double delta_T) {
        return E_modulus * alpha * delta_T;
    }

    double cantilever_beam_deflection(double F, double L, double E_modulus, double I) {
        if (E_modulus <= 0.0 || I <= 0.0) {
            throw std::invalid_argument("cantilever_beam_deflection: modulus and moment of inertia must be positive");
        }
        return (F * L * L * L) / (3.0 * E_modulus * I);
    }

    double continuity_equation_v2(double A1, double v1, double A2) {
        if (A2 <= 0.0) {
            throw std::invalid_argument("continuity_equation_v2: area A2 must be positive");
        }
        return (A1 * v1) / A2;
    }

    double reynolds_number(double rho, double v, double D, double mu) {
        if (mu <= 0.0) {
            throw std::invalid_argument("reynolds_number: dynamic viscosity must be positive");
        }
        return (rho * v * D) / mu;
    }

    double hydrostatic_pressure(double P0, double rho, double g, double h) {
        return P0 + (rho * g * h);
    }

    double buoyant_force(double rho_fluid, double g, double V_displaced) {
        return rho_fluid * g * V_displaced;
    }

    double drag_force(double C_D, double rho, double A, double v) {
        return 0.5 * C_D * rho * A * v * v;
    }

    double lift_force(double C_L, double rho, double A, double v) {
        return 0.5 * C_L * rho * A * v * v;
    }

    double gear_ratio(double N_out, double N_in) {
        if (N_in <= 0.0) {
            throw std::invalid_argument("gear_ratio: input teeth count must be positive");
        }
        return N_out / N_in;
    }

    double newtons_law_cooling_convection(double h_coeff, double A, double delta_T) {
        return h_coeff * A * delta_T;
    }

}