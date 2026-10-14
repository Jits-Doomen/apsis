#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "apsis/ephemerides.hpp"
#include "apsis/propagator.hpp"
#include "apsis/trajectory.hpp"

namespace ApsisPackages::Formulas {
    double projectile_max_heigth(double v0, double theta, double g);
    double projectile_range(double v0, double theta, double g);
    double projectile_time_of_flight(double v0, double theta, double g);
    double average_velocity(double d, double t);
    double average_acceleration(double delta_v, double delta_t);
    double kinematic_final_velocity(double v0, double a, double t);
    double kinematic_position(double x0, double v0, double a, double t);
    double kinematic_velocity_squared(double v0, double a, double delta_x);
    double newton_second_law(double m, double a);
    double weight(double m, double g);
    double gravitational_force(double m1, double m2, double r);
    double linear_momentum(double m, double v);
    double impulse(double F, double delta_t);
    double friction_force(double mu, double N);
    double work_done(double F, double d, double theta);
    double kinetic_energy(double m, double v);
    double gravitational_potential_energy(double m, double g, double h);
    double mechanical_power(double W, double t);
    double hookes_law(double k, double x);
    double elastic_potential_energy(double k, double x);
    double angular_velocity(double delta_theta, double delta_t);
    double angular_acceleration(double delta_omega, double delta_t);
    double tangential_velocity(double omega, double r);
    double centripetal_acceleration(double v, double r);
    double centripetal_force(double m, double v, double r);
    double torque(double r, double F, double theta);
    double moment_of_inertia_point_mass(double m, double r);
    double moment_of_inertia_solid_disk(double m, double r);
    double moment_of_inertia_solid_sphere(double m, double r);
    double moment_of_inertia_thin_rod_center(double m, double L);
    double angular_momentum(double I, double omega);
    double rotational_kinetic_energy(double I, double omega);
    double parallel_axis_theorem(double I_cm, double m, double d);
    double shm_spring_period(double m, double k);
    double simple_pendulum_period(double L, double g);
    double shm_position(double A, double omega, double t, double phi);
    double shm_angular_frequency(double k, double m);
    double damped_oscillation_amplitude(double A0, double b, double m, double t);
    double ideal_gas_pressure(double n, double T, double V);
    double first_law_thermodynamics(double Q, double W);
    double sensible_heat(double m, double c, double delta_T);
    double latent_heat(double m, double L);
    double carnot_efficiency(double T_c, double T_h);
    double entropy_change_reversible(double Q, double T);
    double statistical_entropy(double W);
    double stefan_boltzmann_power(double A, double T);
    double wiens_displacement_law(double T);
    double linear_thermal_expansion(double alpha, double L0, double delta_T);
    double fourier_conduction_rate(double k, double A, double delta_T, double d);
    double rms_molecular_speed(double T, double M_molar);
    double isobaric_work(double P, double delta_V);
    double coulomb_force(double q1, double q2, double r);
    double electric_field_point_charge(double Q_charge, double r);
    double electric_force_on_charge(double q, double E_field);
    double electric_potential_point_charge(double Q_charge, double r);
    double parallel_plate_capacitance(double A, double d);
    double capacitor_stored_energy(double C_cap, double V);
    double ohms_law_voltage(double I, double R);
    double electrical_power(double V, double I);
    double conductor_resistance(double rho, double L, double A);
    double amperes_law_magnetic_field_circulation(double I_enc);
    double faradays_law_emf(double delta_phi_B, double delta_t);
    double self_inductance_emf(double L_ind, double delta_I, double delta_t);
    double inductor_stored_energy(double L_ind, double I);
    double speed_of_em_waves();
    double wave_speed(double f, double lambda);
    double wave_period(double f);
    double snells_law_angle2(double n1, double theta1, double n2);
    double critical_angle(double n1, double n2);
    double doppler_effect_observed_frequency(double f_source, double v_wave, double v_observer, double v_source);
    double thin_lens_image_distance(double f, double d_o);
    double lens_magnification(double d_i, double d_o);
    double wave_intensity(double P, double A);
    double sound_level_decibels(double I);
    double mass_energy_equivalence(double m);
    double lorentz_factor(double v);
    double time_dilation(double delta_t0, double gamma);
    double length_contraction(double L0, double gamma);
    double relativistic_momentum(double m, double v, double gamma);
    double relativistic_total_energy(double m, double gamma);
    double energy_momentum_relation(double p, double m);
    double photon_energy_frequency(double f);
    double de_broglie_wavelength(double p);
    double radioactive_decay_remaining(double N0, double lambda_decay, double t);
    double radioactive_half_life(double lambda_decay);
    double keplers_third_law(double a, double M);
    double stellar_luminosity(double R, double T);
    double flux_luminosity_relation(double L, double d);
    double distance_modulus(double d);
    double parallax_distance(double p_arcsec);
    double hubbles_law_velocity(double H0, double d);
    double cosmological_redshift(double lambda_obs, double lambda_rest);
    double schwarzschild_radius(double M);
    double eddington_luminosity(double M);
    double vis_viva_velocity(double G_const, double M, double r, double a);
    double circular_orbital_velocity(double M, double r);
    double escape_velocity(double M, double r);
    double specific_orbital_energy(double v, double M, double r);
    double orbital_eccentricity_from_apsides(double r_apo, double r_peri);
    double semi_major_axis_from_apsides(double r_apo, double r_peri);
    double hohmann_transfer_delta_v1(double M, double r1, double r2);
    double hohmann_transfer_delta_v2(double M, double r1, double r2);
    double sphere_of_influence_radius(double a, double m_small, double M_large);
    double tsiolkovsky_rocket_equation(double v_e, double m0, double mf);
    double effective_exhaust_velocity(double I_sp, double g0_acc);
    double specific_impulse(double F_thrust, double m_dot, double g0_acc);
    double thrust_equation(double m_dot, double v_e, double p_e, double p_a, double A_e);
    double rocket_mass_ratio(double m0, double mf);
    double propellant_mass_fraction(double m_p, double m0);
    double thrust_to_weight_ratio(double F_thrust, double m, double g0_acc);
    double rocket_burn_time(double m_p, double m_dot);
    double normal_stress(double F, double A);
    double normal_strain(double delta_L, double L0);
    double youngs_modulus(double sigma, double epsilon);
    double shear_modulus(double tau, double gamma_strain);
    double poissons_ratio(double epsilon_transverse, double epsilon_axial);
    double bending_stress(double M_moment, double y, double I);
    double torsional_shear_stress(double T_torque, double r, double J);
    double angle_of_twist(double T_torque, double L, double J, double G_modulus);
    double euler_buckling_load(double E_modulus, double I, double K_factor, double L);
    double factor_of_safety(double sigma_yield, double sigma_allowable);
    double thermal_stress_constrained(double E_modulus, double alpha, double delta_T);
    double cantilever_beam_deflection(double F, double L, double E_modulus, double I);
    double continuity_equation_v2(double A1, double v1, double A2);
    double reynolds_number(double rho, double v, double D, double mu);
    double hydrostatic_pressure(double P0, double rho, double g, double h);
    double buoyant_force(double rho_fluid, double g, double V_displaced);
    double drag_force(double C_D, double rho, double A, double v);
    double lift_force(double C_L, double rho, double A, double v);
    double gear_ratio(double N_out, double N_in);
    double newtons_law_cooling_convection(double h_coeff, double A, double delta_T);
} // namespace ApsisPackages::Formulas

// Apsis: a small C++20 interpreter for scientific and astrodynamics scripts.
// Build: CMake / C++20. Run: apsis path/to/script.jdx

namespace apsis {

    struct ApsisError : std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    struct Dimension {
        // SI base dimensions: length, mass, time, electric current, temperature
        std::array<int, 5> e{ 0, 0, 0, 0, 0 };

        bool operator==(const Dimension& other) const { return e == other.e; }
        bool dimensionless() const {
            return std::all_of(e.begin(), e.end(), [](int x) { return x == 0; });
        }
        Dimension operator+(const Dimension& b) const {
            Dimension r;
            for (std::size_t i = 0; i < e.size(); ++i) r.e[i] = e[i] + b.e[i];
            return r;
        }
        Dimension operator-(const Dimension& b) const {
            Dimension r;
            for (std::size_t i = 0; i < e.size(); ++i) r.e[i] = e[i] - b.e[i];
            return r;
        }
        Dimension scaled(int n) const {
            Dimension r;
            for (std::size_t i = 0; i < e.size(); ++i) r.e[i] = e[i] * n;
            return r;
        }
    };

    static std::string dimensionName(const Dimension& d) {
        static const std::array<const char*, 5> names{ "m", "kg", "s", "A", "K" };
        if (d.dimensionless()) return "1";

        std::ostringstream numerator;
        std::ostringstream denominator;
        bool hasNumerator = false;
        bool hasDenominator = false;

        for (std::size_t i = 0; i < d.e.size(); ++i) {
            const int exponent = d.e[i];
            if (exponent == 0) continue;

            std::ostringstream& out = exponent > 0 ? numerator : denominator;
            bool& hasTerm = exponent > 0 ? hasNumerator : hasDenominator;

            if (hasTerm) out << '*';
            hasTerm = true;
            out << names[i];

            const int magnitude = exponent > 0 ? exponent : -exponent;
            if (magnitude != 1) out << '^' << magnitude;
        }

        std::string result = hasNumerator ? numerator.str() : "1";
        if (hasDenominator) result += "/" + denominator.str();
        return result;
    }

    struct Quantity {
        double value = 0.0;
        Dimension dim{};

        Quantity() = default;
        Quantity(double v, Dimension d = {}) : value(v), dim(d) {}

        Quantity add(const Quantity& b) const {
            if (!(dim == b.dim)) throw ApsisError("dimension mismatch: cannot add " +
                dimensionName(dim) + " and " + dimensionName(b.dim));
            return { value + b.value, dim };
        }
        Quantity sub(const Quantity& b) const {
            if (!(dim == b.dim)) throw ApsisError("dimension mismatch: cannot subtract " +
                dimensionName(dim) + " and " + dimensionName(b.dim));
            return { value - b.value, dim };
        }
        Quantity mul(const Quantity& b) const { return { value * b.value, dim + b.dim }; }
        Quantity div(const Quantity& b) const {
            if (b.value == 0.0) throw ApsisError("division by zero");
            return { value / b.value, dim - b.dim };
        }
    };

    struct Vec3 {
        Quantity x, y, z;

        Vec3(Quantity a, Quantity b, Quantity c) : x(a), y(b), z(c) {
            if (!(x.dim == y.dim && y.dim == z.dim))
                throw ApsisError("all components of a vector must have matching dimensions");
        }
        Quantity norm() const {
            return { std::sqrt(x.value * x.value + y.value * y.value + z.value * z.value), x.dim };
        }
        Vec3 unit() const {
            const double n = norm().value;
            if (n == 0.0) throw ApsisError("cannot normalize a zero vector");
            return { {x.value / n, {}}, {y.value / n, {}}, {z.value / n, {}} };
        }
        Quantity dot(const Vec3& b) const {
            return { x.value * b.x.value + y.value * b.y.value + z.value * b.z.value, x.dim + b.x.dim };
        }
        Vec3 cross(const Vec3& b) const {
            return {
                {y.value * b.z.value - z.value * b.y.value, y.dim + b.z.dim},
                {z.value * b.x.value - x.value * b.z.value, z.dim + b.x.dim},
                {x.value * b.y.value - y.value * b.x.value, x.dim + b.x.dim}
            };
        }
    };

    struct Function;
    struct Environment;
    using FunctionPtr = std::shared_ptr<Function>;
    using Value = std::variant<std::monostate, bool, std::string, Quantity, Vec3, FunctionPtr>;

    static std::string numberText(double v) {
        std::ostringstream out;
        out << std::setprecision(12) << v;
        return out.str();
    }
    static std::string valueText(const Value& v) {
        if (std::holds_alternative<std::monostate>(v)) return "nothing";
        if (auto p = std::get_if<bool>(&v)) return *p ? "true" : "false";
        if (auto p = std::get_if<std::string>(&v)) return *p;
        if (auto p = std::get_if<Quantity>(&v)) {
            return numberText(p->value) + (p->dim.dimensionless() ? "" : " " + dimensionName(p->dim));
        }
        if (auto p = std::get_if<Vec3>(&v)) {
            return "vec(" + valueText(Value{ p->x }) + ", " + valueText(Value{ p->y }) + ", " +
                valueText(Value{ p->z }) + ")";
        }
        return "<function>";
    }
    static Quantity asQuantity(const Value& v, const std::string& where) {
        if (auto p = std::get_if<Quantity>(&v)) return *p;
        throw ApsisError(where + " requires a quantity");
    }
    static bool truthy(const Value& v) {
        if (auto p = std::get_if<bool>(&v)) return *p;
        if (auto p = std::get_if<Quantity>(&v)) {
            if (!p->dim.dimensionless()) throw ApsisError("condition must be dimensionless");
            return p->value != 0.0;
        }
        if (std::holds_alternative<std::monostate>(v)) return false;
        throw ApsisError("condition must be a boolean or dimensionless number");
    }

    enum class TokenType {
        End, Number, String, Identifier,
        LParen, RParen, LBrace, RBrace, Comma, Semicolon, Newline,
        Plus, Minus, Star, Slash, Caret, Percent,
        Bang, Equal, EqualEqual, BangEqual, Less, LessEqual, Greater, GreaterEqual,
        AndAnd, OrOr,
        Define, Const, Fn, Return, If, Else, While, Import, True, False
    };

    struct Token {
        TokenType type;
        std::string text;
        double number = 0.0;
        int line = 1;
        int column = 1;
    };

    class Lexer {
        const std::string& source;
        std::size_t pos = 0;
        int line = 1, col = 1;
        std::vector<Token> tokens;

        char peek(std::size_t n = 0) const {
            return pos + n < source.size() ? source[pos + n] : '\0';
        }
        char advance() {
            char c = peek();
            if (c == '\0') return c;
            ++pos;
            if (c == '\n') { ++line; col = 1; }
            else ++col;
            return c;
        }
        void add(TokenType t, std::string s, int l, int c, double n = 0.0) {
            tokens.push_back(Token{ t, std::move(s), n, l, c });
        }
        [[noreturn]] void fail(const std::string& message, int l, int c) const {
            throw ApsisError("line " + std::to_string(l) + ", column " +
                std::to_string(c) + ": " + message);
        }
    public:
        explicit Lexer(const std::string& src) : source(src) {}

        std::vector<Token> scan() {
            while (peek() != '\0') {
                char c = peek();
                if (c == ' ' || c == '\t' || c == '\r') { advance(); continue; }
                int l = line, cc = col;
                if (c == '\n') { advance(); add(TokenType::Newline, "\\n", l, cc); continue; }
                if (c == '#' || (c == '/' && peek(1) == '/')) {
                    while (peek() != '\0' && peek() != '\n') advance();
                    continue;
                }
                if (std::isdigit(static_cast<unsigned char>(c)) ||
                    (c == '.' && std::isdigit(static_cast<unsigned char>(peek(1))))) {
                    std::size_t start = pos;
                    bool seenDot = false;
                    if (c == '.') { seenDot = true; advance(); }
                    while (std::isdigit(static_cast<unsigned char>(peek()))) advance();
                    if (!seenDot && peek() == '.') {
                        advance();
                        while (std::isdigit(static_cast<unsigned char>(peek()))) advance();
                    }
                    if (peek() == 'e' || peek() == 'E') {
                        advance();
                        if (peek() == '+' || peek() == '-') advance();
                        if (!std::isdigit(static_cast<unsigned char>(peek())))
                            fail("malformed scientific-notation number", l, cc);
                        while (std::isdigit(static_cast<unsigned char>(peek()))) advance();
                    }
                    std::string s = source.substr(start, pos - start);
                    char* end = nullptr;
                    double n = std::strtod(s.c_str(), &end);
                    if (!end || *end != '\0' || !std::isfinite(n))
                        fail("invalid numeric literal", l, cc);
                    add(TokenType::Number, s, l, cc, n);
                    continue;
                }
                if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
                    std::size_t start = pos;
                    while (std::isalnum(static_cast<unsigned char>(peek())) || peek() == '_') advance();
                    std::string s = source.substr(start, pos - start);
                    static const std::unordered_map<std::string, TokenType> keywords{
                        {"define",TokenType::Define},{"const",TokenType::Const},{"fn",TokenType::Fn},
                        {"return",TokenType::Return},{"if",TokenType::If},{"else",TokenType::Else},
                        {"while",TokenType::While},{"import",TokenType::Import},
                        {"true",TokenType::True},{"false",TokenType::False}
                    };
                    auto it = keywords.find(s);
                    add(it == keywords.end() ? TokenType::Identifier : it->second, s, l, cc);
                    continue;
                }
                if (c == '"') {
                    advance();
                    std::string s;
                    while (peek() != '"' && peek() != '\0') {
                        char ch = advance();
                        if (ch == '\n') fail("newline in string literal", l, cc);
                        if (ch == '\\') {
                            char e = advance();
                            switch (e) {
                            case 'n': s.push_back('\n'); break;
                            case 't': s.push_back('\t'); break;
                            case '"': s.push_back('"'); break;
                            case '\\': s.push_back('\\'); break;
                            default: fail("unknown string escape", line, col);
                            }
                        }
                        else s.push_back(ch);
                    }
                    if (peek() != '"') fail("unterminated string literal", l, cc);
                    advance();
                    add(TokenType::String, s, l, cc);
                    continue;
                }
                auto one = [&](TokenType t) { advance(); add(t, std::string(1, c), l, cc); };
                auto two = [&](char second, TokenType yes, TokenType no) {
                    advance();
                    if (peek() == second) { advance(); add(yes, std::string{ c,second }, l, cc); }
                    else add(no, std::string(1, c), l, cc);
                    };
                switch (c) {
                case '(': one(TokenType::LParen); break;
                case ')': one(TokenType::RParen); break;
                case '{': one(TokenType::LBrace); break;
                case '}': one(TokenType::RBrace); break;
                case ',': one(TokenType::Comma); break;
                case ';': one(TokenType::Semicolon); break;
                case '+': one(TokenType::Plus); break;
                case '-': one(TokenType::Minus); break;
                case '*': one(TokenType::Star); break;
                case '/': one(TokenType::Slash); break;
                case '^': one(TokenType::Caret); break;
                case '%': one(TokenType::Percent); break;
                case '!': two('=', TokenType::BangEqual, TokenType::Bang); break;
                case '=': two('=', TokenType::EqualEqual, TokenType::Equal); break;
                case '<': two('=', TokenType::LessEqual, TokenType::Less); break;
                case '>': two('=', TokenType::GreaterEqual, TokenType::Greater); break;
                case '&':
                    advance();
                    if (peek() != '&') fail("expected '&' after '&'", l, cc);
                    advance(); add(TokenType::AndAnd, "&&", l, cc); break;
                case '|':
                    advance();
                    if (peek() != '|') fail("expected '|' after '|'", l, cc);
                    advance(); add(TokenType::OrOr, "||", l, cc); break;
                default: fail("unexpected character '" + std::string(1, c) + "'", l, cc);
                }
            }
            tokens.push_back(Token{ TokenType::End, "", 0.0, line, col });
            return tokens;
        }
    };

    struct Expr {
        virtual ~Expr() = default;
        virtual Value eval(Environment&) = 0;
    };
    struct Stmt {
        virtual ~Stmt() = default;
        virtual void exec(Environment&) = 0;
    };
    using ExprPtr = std::unique_ptr<Expr>;
    using StmtPtr = std::unique_ptr<Stmt>;

    struct ReturnSignal { Value value; };

    static void installPackage(Environment& env, const std::string& packageName);

    struct Environment : std::enable_shared_from_this<Environment> {
        std::unordered_map<std::string, Value> values;
        std::unordered_map<std::string, bool> immutable;
        std::shared_ptr<Environment> parent;

        explicit Environment(std::shared_ptr<Environment> p = nullptr) : parent(std::move(p)) {}
        void define(const std::string& name, Value value, bool constant = false) {
            if (values.contains(name)) throw ApsisError("name '" + name + "' is already defined in this scope");
            values[name] = std::move(value);
            immutable[name] = constant;
        }
        std::pair<Environment*, Value*> locate(const std::string& name) {
            auto it = values.find(name);
            if (it != values.end()) return { this, &it->second };
            if (parent) return parent->locate(name);
            throw ApsisError("unknown name '" + name + "'");
        }
        Value get(const std::string& name) { return *locate(name).second; }
        void assign(const std::string& name, Value value) {
            auto [env, slot] = locate(name);
            if (env->immutable[name]) throw ApsisError("cannot assign to constant '" + name + "'");
            *slot = std::move(value);
        }
    };

    struct LiteralExpr : Expr {
        Value value;
        explicit LiteralExpr(Value v) : value(std::move(v)) {}
        Value eval(Environment&) override { return value; }
    };
    struct NameExpr : Expr {
        std::string name;
        explicit NameExpr(std::string n) : name(std::move(n)) {}
        Value eval(Environment& e) override { return e.get(name); }
    };
    struct AssignExpr : Expr {
        std::string name; ExprPtr value;
        AssignExpr(std::string n, ExprPtr v) : name(std::move(n)), value(std::move(v)) {}
        Value eval(Environment& e) override {
            Value v = value->eval(e); e.assign(name, v); return v;
        }
    };
    struct UnaryExpr : Expr {
        TokenType op; ExprPtr right;
        UnaryExpr(TokenType o, ExprPtr r) : op(o), right(std::move(r)) {}
        Value eval(Environment& e) override {
            Value v = right->eval(e);
            if (op == TokenType::Bang) return Value{ !truthy(v) };
            if (op == TokenType::Minus) {
                if (auto q = std::get_if<Quantity>(&v)) return Quantity{ -q->value, q->dim };
                if (auto p = std::get_if<Vec3>(&v))
                    return Vec3{ {-p->x.value,p->x.dim},{-p->y.value,p->y.dim},{-p->z.value,p->z.dim} };
                throw ApsisError("unary '-' requires a quantity or vector");
            }
            if (op == TokenType::Plus) {
                if (std::holds_alternative<Quantity>(v)) return v;
                throw ApsisError("unary '+' requires a quantity");
            }
            throw ApsisError("unsupported unary operator");
        }
    };

    static bool equalValues(const Value& a, const Value& b) {
        if (a.index() != b.index()) return false;
        if (auto x = std::get_if<bool>(&a)) return *x == std::get<bool>(b);
        if (auto x = std::get_if<std::string>(&a)) return *x == std::get<std::string>(b);
        if (std::holds_alternative<std::monostate>(a)) return true;
        if (auto x = std::get_if<Quantity>(&a)) {
            auto y = std::get<Quantity>(b);
            if (!(x->dim == y.dim)) throw ApsisError("cannot compare quantities with different dimensions");
            return x->value == y.value;
        }
        if (auto x = std::get_if<Vec3>(&a)) {
            auto y = std::get<Vec3>(b);
            return x->x.value == y.x.value && x->y.value == y.y.value && x->z.value == y.z.value &&
                x->x.dim == y.x.dim;
        }
        return false;
    }

    struct BinaryExpr : Expr {
        ExprPtr left; TokenType op; ExprPtr right;
        BinaryExpr(ExprPtr a, TokenType o, ExprPtr b) : left(std::move(a)), op(o), right(std::move(b)) {}
        Value eval(Environment& e) override {
            Value a = left->eval(e);
            if (op == TokenType::AndAnd) return Value{ truthy(a) && truthy(right->eval(e)) };
            if (op == TokenType::OrOr) return Value{ truthy(a) || truthy(right->eval(e)) };
            Value b = right->eval(e);
            if (op == TokenType::EqualEqual) return Value{ equalValues(a,b) };
            if (op == TokenType::BangEqual) return Value{ !equalValues(a,b) };
            if (auto av = std::get_if<Quantity>(&a)) {
                auto bv = std::get_if<Quantity>(&b);
                if (!bv) throw ApsisError("operator requires two quantities");
                switch (op) {
                case TokenType::Plus: return av->add(*bv);
                case TokenType::Minus: return av->sub(*bv);
                case TokenType::Star: return av->mul(*bv);
                case TokenType::Slash: return av->div(*bv);
                case TokenType::Percent:
                    if (!(av->dim == bv->dim)) throw ApsisError("remainder requires matching dimensions");
                    if (bv->value == 0.0) throw ApsisError("remainder by zero");
                    return Quantity{ std::fmod(av->value,bv->value),av->dim };
                case TokenType::Caret: {
                    if (!bv->dim.dimensionless()) throw ApsisError("exponent must be dimensionless");
                    double rounded = std::round(bv->value);
                    if (!av->dim.dimensionless() && std::abs(bv->value - rounded) > 1e-12)
                        throw ApsisError("dimensioned quantities require an integer exponent");
                    if (av->value < 0 && std::floor(bv->value) != bv->value)
                        throw ApsisError("fractional power of a negative value");
                    if (av->value == 0 && bv->value < 0) throw ApsisError("zero cannot have a negative exponent");
                    if (std::abs(rounded) > 1000 && !av->dim.dimensionless())
                        throw ApsisError("exponent is too large for dimensional arithmetic");
                    Dimension d = av->dim;
                    if (!d.dimensionless()) d = d.scaled(static_cast<int>(rounded));
                    return Quantity{ std::pow(av->value,bv->value),d };
                }
                case TokenType::Less: case TokenType::LessEqual:
                case TokenType::Greater: case TokenType::GreaterEqual:
                    if (!(av->dim == bv->dim)) throw ApsisError("comparison requires matching dimensions");
                    if (op == TokenType::Less) return Value{ av->value < bv->value };
                    if (op == TokenType::LessEqual) return Value{ av->value <= bv->value };
                    if (op == TokenType::Greater) return Value{ av->value > bv->value };
                    return Value{ av->value >= bv->value };
                default: break;
                }
            }
            if (auto va = std::get_if<Vec3>(&a)) {
                if (auto vb = std::get_if<Vec3>(&b)) {
                    if (op == TokenType::Plus || op == TokenType::Minus) {
                        Quantity x = op == TokenType::Plus ? va->x.add(vb->x) : va->x.sub(vb->x);
                        Quantity y = op == TokenType::Plus ? va->y.add(vb->y) : va->y.sub(vb->y);
                        Quantity z = op == TokenType::Plus ? va->z.add(vb->z) : va->z.sub(vb->z);
                        return Vec3{ x,y,z };
                    }
                    throw ApsisError("unsupported vector operator");
                }
                if (auto qb = std::get_if<Quantity>(&b)) {
                    if (op == TokenType::Star) return Vec3{ va->x.mul(*qb),va->y.mul(*qb),va->z.mul(*qb) };
                    if (op == TokenType::Slash) return Vec3{ va->x.div(*qb),va->y.div(*qb),va->z.div(*qb) };
                }
            }
            if (auto qa = std::get_if<Quantity>(&a)) {
                if (auto vb = std::get_if<Vec3>(&b); vb && op == TokenType::Star)
                    return Vec3{ vb->x.mul(*qa),vb->y.mul(*qa),vb->z.mul(*qa) };
            }
            if (op == TokenType::Plus) {
                if (auto sa = std::get_if<std::string>(&a)) {
                    if (auto sb = std::get_if<std::string>(&b)) return *sa + *sb;
                }
            }
            throw ApsisError("invalid operands for operator '" + std::string(
                op == TokenType::Plus ? "+" : op == TokenType::Minus ? "-" :
                op == TokenType::Star ? "*" : op == TokenType::Slash ? "/" :
                op == TokenType::Caret ? "^" : "?") + "'");
        }
    };

    struct Function {
        std::vector<std::string> parameters;
        std::vector<StmtPtr> body;
        std::shared_ptr<Environment> closure;
        bool builtin = false;
        std::function<Value(const std::vector<Value>&)> native;
    };

    struct CallExpr : Expr {
        ExprPtr callee; std::vector<ExprPtr> args;
        CallExpr(ExprPtr c, std::vector<ExprPtr> a) : callee(std::move(c)), args(std::move(a)) {}
        Value eval(Environment& e) override;
    };

    struct ImportStmt : Stmt {
        std::string packageName;
        explicit ImportStmt(std::string name) : packageName(std::move(name)) {}
        void exec(Environment& e) override { installPackage(e, packageName); }
    };

    struct ExprStmt : Stmt {
        ExprPtr expr;
        explicit ExprStmt(ExprPtr x) : expr(std::move(x)) {}
        void exec(Environment& e) override { (void)expr->eval(e); }
    };
    struct DefineStmt : Stmt {
        std::string name; ExprPtr initializer; bool constant;
        DefineStmt(std::string n, ExprPtr x, bool c) : name(std::move(n)), initializer(std::move(x)), constant(c) {}
        void exec(Environment& e) override { e.define(name, initializer->eval(e), constant); }
    };
    struct BlockStmt : Stmt {
        std::vector<StmtPtr> statements;
        explicit BlockStmt(std::vector<StmtPtr> s) : statements(std::move(s)) {}
        void exec(Environment& e) override {
            // Blocks execute in the current environment. Function calls establish their own scope.
            // Keeping block execution in-place also allows nested blocks without dangling scopes.
            for (auto& s : statements) s->exec(e);
        }
    };
    struct IfStmt : Stmt {
        ExprPtr condition; StmtPtr yes, no;
        IfStmt(ExprPtr c, StmtPtr y, StmtPtr n) : condition(std::move(c)), yes(std::move(y)), no(std::move(n)) {}
        void exec(Environment& e) override {
            if (truthy(condition->eval(e))) yes->exec(e);
            else if (no) no->exec(e);
        }
    };
    struct WhileStmt : Stmt {
        ExprPtr condition; StmtPtr body;
        WhileStmt(ExprPtr c, StmtPtr b) : condition(std::move(c)), body(std::move(b)) {}
        void exec(Environment& e) override {
            std::size_t iterations = 0;
            while (truthy(condition->eval(e))) {
                if (++iterations > 10'000'000) throw ApsisError("loop exceeded 10,000,000 iterations");
                body->exec(e);
            }
        }
    };
    struct FunctionStmt : Stmt {
        std::string name; std::vector<std::string> params; std::vector<StmtPtr> body;
        FunctionStmt(std::string n, std::vector<std::string> p, std::vector<StmtPtr> b)
            : name(std::move(n)), params(std::move(p)), body(std::move(b)) {
        }
        void exec(Environment& e) override {
            auto f = std::make_shared<Function>();
            f->parameters = params;
            f->body = std::move(body);
            f->closure = e.shared_from_this();
            e.define(name, f, true);
        }
    };
    struct ReturnStmt : Stmt {
        ExprPtr value;
        explicit ReturnStmt(ExprPtr v) : value(std::move(v)) {}
        void exec(Environment& e) override { throw ReturnSignal{ value ? value->eval(e) : Value{} }; }
    };
    struct PrintStmt : Stmt {
        std::vector<ExprPtr> values;
        explicit PrintStmt(std::vector<ExprPtr> v) : values(std::move(v)) {}
        void exec(Environment& e) override {
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (i) std::cout << ' ';
                std::cout << valueText(values[i]->eval(e));
            }
            std::cout << '\n';
        }
    };

    Value CallExpr::eval(Environment& e) {
        Value cv = callee->eval(e);
        auto fp = std::get_if<FunctionPtr>(&cv);
        if (!fp || !*fp) throw ApsisError("attempted to call a non-function value");
        std::vector<Value> argv;
        for (auto& a : args) argv.push_back(a->eval(e));
        Function& f = **fp;
        if (f.builtin) return f.native(argv);
        if (argv.size() != f.parameters.size())
            throw ApsisError("function expects " + std::to_string(f.parameters.size()) +
                " argument(s), received " + std::to_string(argv.size()));
        auto local = std::make_shared<Environment>(f.closure);
        for (std::size_t i = 0; i < argv.size(); ++i) local->define(f.parameters[i], argv[i]);
        try {
            for (auto& s : f.body) s->exec(*local);
        }
        catch (const ReturnSignal& r) { return r.value; }
        return {};
    }

    class Parser {
        std::vector<Token> ts;
        std::size_t current = 0;

        const Token& peek() const { return ts[current]; }
        const Token& previous() const { return ts[current - 1]; }
        bool atEnd() const { return peek().type == TokenType::End; }
        bool check(TokenType t) const { return peek().type == t; }
        Token advance() { if (!atEnd()) ++current; return previous(); }
        bool match(std::initializer_list<TokenType> types) {
            for (auto t : types) if (check(t)) { advance(); return true; }
            return false;
        }
        Token consume(TokenType t, const std::string& msg) {
            if (check(t)) return advance();
            fail(msg);
        }
        [[noreturn]] void fail(const std::string& msg) const {
            const auto& t = peek();
            throw ApsisError("line " + std::to_string(t.line) + ", column " +
                std::to_string(t.column) + ": " + msg);
        }
        void separators() { while (match({ TokenType::Newline,TokenType::Semicolon })) {} }
        void endStatement() {
            if (!check(TokenType::RBrace) && !check(TokenType::End) &&
                !check(TokenType::Newline) && !check(TokenType::Semicolon))
                fail("expected end of statement");
            separators();
        }

        StmtPtr statement() {
            if (match({ TokenType::Newline,TokenType::Semicolon })) {
                while (match({ TokenType::Newline,TokenType::Semicolon })) {}
                return std::make_unique<ExprStmt>(std::make_unique<LiteralExpr>(Value{}));
            }
            if (match({ TokenType::Import })) {
                Token package = consume(TokenType::Identifier, "expected package name after 'import'");
                endStatement();
                return std::make_unique<ImportStmt>(package.text);
            }
            if (match({ TokenType::Define,TokenType::Const })) {
                bool isConst = previous().type == TokenType::Const;
                Token name = consume(TokenType::Identifier, "expected a name after declaration keyword");
                consume(TokenType::Equal, "expected '=' after declaration name");
                auto value = expression();
                endStatement();
                return std::make_unique<DefineStmt>(name.text, std::move(value), isConst);
            }
            if (match({ TokenType::Fn })) {
                Token name = consume(TokenType::Identifier, "expected function name");
                consume(TokenType::LParen, "expected '(' after function name");
                std::vector<std::string> params;
                if (!check(TokenType::RParen)) {
                    do {
                        if (params.size() >= 255) fail("functions may have at most 255 parameters");
                        params.push_back(consume(TokenType::Identifier, "expected parameter name").text);
                    } while (match({ TokenType::Comma }));
                }
                consume(TokenType::RParen, "expected ')' after parameters");
                separators();
                consume(TokenType::LBrace, "expected '{' before function body");
                auto body = blockContents();
                return std::make_unique<FunctionStmt>(name.text, std::move(params), std::move(body));
            }
            if (match({ TokenType::Return })) {
                ExprPtr v;
                if (!check(TokenType::Newline) && !check(TokenType::Semicolon) &&
                    !check(TokenType::RBrace) && !check(TokenType::End)) v = expression();
                endStatement();
                return std::make_unique<ReturnStmt>(std::move(v));
            }
            if (match({ TokenType::If })) {
                consume(TokenType::LParen, "expected '(' after if");
                auto condition = expression();
                consume(TokenType::RParen, "expected ')' after condition");
                separators();
                auto yes = statementOrBlock();
                separators();
                StmtPtr no;
                if (match({ TokenType::Else })) {
                    separators();
                    no = statementOrBlock();
                }
                return std::make_unique<IfStmt>(std::move(condition), std::move(yes), std::move(no));
            }
            if (match({ TokenType::While })) {
                consume(TokenType::LParen, "expected '(' after while");
                auto condition = expression();
                consume(TokenType::RParen, "expected ')' after condition");
                separators();
                return std::make_unique<WhileStmt>(std::move(condition), statementOrBlock());
            }
            if (match({ TokenType::LBrace })) return std::make_unique<BlockStmt>(blockContents());
            if (check(TokenType::Identifier) && peek().text == "print") {
                advance();
                consume(TokenType::LParen, "expected '(' after print");
                std::vector<ExprPtr> args;
                if (!check(TokenType::RParen)) {
                    do { args.push_back(expression()); } while (match({ TokenType::Comma }));
                }
                consume(TokenType::RParen, "expected ')' after print arguments");
                endStatement();
                return std::make_unique<PrintStmt>(std::move(args));
            }
            auto e = expression();
            endStatement();
            return std::make_unique<ExprStmt>(std::move(e));
        }
        StmtPtr statementOrBlock() {
            if (match({ TokenType::LBrace })) return std::make_unique<BlockStmt>(blockContents());
            return statement();
        }
        std::vector<StmtPtr> blockContents() {
            std::vector<StmtPtr> body;
            separators();
            while (!check(TokenType::RBrace) && !atEnd()) {
                body.push_back(statement());
                separators();
            }
            consume(TokenType::RBrace, "expected '}' after block");
            return body;
        }
        ExprPtr expression() { return assignment(); }
        ExprPtr assignment() {
            auto expr = logicalOr();
            if (match({ TokenType::Equal })) {
                auto rhs = assignment();
                if (auto n = dynamic_cast<NameExpr*>(expr.get())) {
                    std::string name = n->name;
                    return std::make_unique<AssignExpr>(name, std::move(rhs));
                }
                fail("assignment target must be a variable name");
            }
            return expr;
        }
        ExprPtr logicalOr() {
            auto e = logicalAnd();
            while (match({ TokenType::OrOr })) { auto op = previous().type; auto rhs = logicalAnd(); e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs)); }
            return e;
        }
        ExprPtr logicalAnd() {
            auto e = equality();
            while (match({ TokenType::AndAnd })) { auto op = previous().type; auto rhs = equality(); e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs)); }
            return e;
        }
        ExprPtr equality() {
            auto e = comparison();
            while (match({ TokenType::EqualEqual,TokenType::BangEqual })) {
                auto op = previous().type; auto rhs = comparison();
                e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs));
            }
            return e;
        }
        ExprPtr comparison() {
            auto e = term();
            while (match({ TokenType::Less,TokenType::LessEqual,TokenType::Greater,TokenType::GreaterEqual })) {
                auto op = previous().type; auto rhs = term();
                e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs));
            }
            return e;
        }
        ExprPtr term() {
            auto e = factor();
            while (match({ TokenType::Plus,TokenType::Minus })) {
                auto op = previous().type; auto rhs = factor();
                e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs));
            }
            return e;
        }
        ExprPtr factor() {
            auto e = unary();
            while (match({ TokenType::Star,TokenType::Slash,TokenType::Percent })) {
                auto op = previous().type; auto rhs = unary();
                e = std::make_unique<BinaryExpr>(std::move(e), op, std::move(rhs));
            }
            return e;
        }
        ExprPtr unary() {
            if (match({ TokenType::Bang,TokenType::Minus,TokenType::Plus }))
                return std::make_unique<UnaryExpr>(previous().type, unary());
            return power();
        }
        ExprPtr power() {
            auto e = call();
            if (match({ TokenType::Caret })) e = std::make_unique<BinaryExpr>(std::move(e), TokenType::Caret, unary());
            return e;
        }
        ExprPtr call() {
            auto e = primary();
            while (match({ TokenType::LParen })) {
                std::vector<ExprPtr> args;
                if (!check(TokenType::RParen)) {
                    do {
                        if (args.size() >= 255) fail("calls may have at most 255 arguments");
                        args.push_back(expression());
                    } while (match({ TokenType::Comma }));
                }
                consume(TokenType::RParen, "expected ')' after arguments");
                e = std::make_unique<CallExpr>(std::move(e), std::move(args));
            }
            return e;
        }
        ExprPtr primary() {
            if (match({ TokenType::Number })) return std::make_unique<LiteralExpr>(Quantity{ previous().number });
            if (match({ TokenType::String })) return std::make_unique<LiteralExpr>(previous().text);
            if (match({ TokenType::True })) return std::make_unique<LiteralExpr>(true);
            if (match({ TokenType::False })) return std::make_unique<LiteralExpr>(false);
            if (match({ TokenType::Identifier })) return std::make_unique<NameExpr>(previous().text);
            if (match({ TokenType::LParen })) {
                auto e = expression();
                consume(TokenType::RParen, "expected ')' after expression");
                return e;
            }
            fail("expected an expression");
        }
    public:
        explicit Parser(std::vector<Token> tokens) : ts(std::move(tokens)) {}
        std::vector<StmtPtr> parse() {
            std::vector<StmtPtr> result;
            separators();
            while (!atEnd()) {
                result.push_back(statement());
                separators();
            }
            return result;
        }
    };

    static void requireCount(const std::vector<Value>& a, std::size_t n, const std::string& name) {
        if (a.size() != n) throw ApsisError(name + " expects " + std::to_string(n) + " argument(s)");
    }
    static Quantity qarg(const std::vector<Value>& a, std::size_t i, const std::string& name) {
        return asQuantity(a.at(i), name);
    }
    static void requireDimensionless(const Quantity& q, const std::string& name) {
        if (!q.dim.dimensionless()) throw ApsisError(name + " requires a dimensionless argument");
    }
    static Value nativeFunction(std::function<Value(const std::vector<Value>&)> fn) {
        auto f = std::make_shared<Function>();
        f->builtin = true;
        f->native = std::move(fn);
        return f;
    }
    static void installBuiltins(const std::shared_ptr<Environment>& env) {
        auto addFn = [&](const std::string& n, auto f) { env->define(n, nativeFunction(f), true); };
        auto addQ = [&](const std::string& n, double v, Dimension d = {}) { env->define(n, Quantity{ v,d }, true); };
        auto dim = [](int l, int m, int t, int i = 0, int k = 0) {
            Dimension d; d.e = { l,m,t,i,k }; return d;
            };

        addQ("PI", std::acos(-1.0));
        addQ("G", 6.67430e-11, dim(3, -1, -2));
        addQ("c", 299792458.0, dim(1, 0, -1));
        addQ("AU", 149597870700.0, dim(1, 0, 0));
        addQ("M_sun", 1.98847e30, dim(0, 1, 0));
        addQ("mu_sun", 1.32712440018e20, dim(3, 0, -2));
        addQ("M_earth", 5.9722e24, dim(0, 1, 0));
        addQ("R_earth", 6371000.0, dim(1, 0, 0));
        addQ("mu_earth", 3.986004418e14, dim(3, 0, -2));

        addQ("m", 1.0, dim(1, 0, 0)); addQ("km", 1000.0, dim(1, 0, 0));
        addQ("cm", 0.01, dim(1, 0, 0)); addQ("mm", 0.001, dim(1, 0, 0));
        addQ("kg", 1.0, dim(0, 1, 0)); addQ("g", 0.001, dim(0, 1, 0));
        addQ("s", 1.0, dim(0, 0, 1)); addQ("min", 60.0, dim(0, 0, 1));
        addQ("hr", 3600.0, dim(0, 0, 1)); addQ("day", 86400.0, dim(0, 0, 1));
        addQ("A", 1.0, dim(0, 0, 0, 1)); addQ("K", 1.0, dim(0, 0, 0, 0, 1));
        addQ("N", 1.0, dim(1, 1, -2)); addQ("J", 1.0, dim(2, 1, -2));
        addQ("W", 1.0, dim(2, 1, -3)); addQ("Pa", 1.0, dim(-1, 1, -2));
        addQ("Hz", 1.0, dim(0, 0, -1));

        addFn("sqrt", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "sqrt"); auto q = qarg(a, 0, "sqrt");
            for (int e : q.dim.e) if (e % 2 != 0) throw ApsisError("sqrt requires even dimension exponents");
            if (q.value < 0) throw ApsisError("sqrt requires a non-negative value");
            Dimension d; for (std::size_t i = 0;i < 5;++i) d.e[i] = q.dim.e[i] / 2;
            return Quantity{ std::sqrt(q.value),d };
            });
        auto trig = [&](const std::string& name, auto fn) {
            addFn(name, [name, fn](const std::vector<Value>& a)->Value {
                requireCount(a, 1, name); auto q = qarg(a, 0, name); requireDimensionless(q, name);
                return Quantity{ fn(q.value) };
                });
            };
        trig("sin", static_cast<double(*)(double)>(std::sin));
        trig("cos", static_cast<double(*)(double)>(std::cos));
        trig("tan", static_cast<double(*)(double)>(std::tan));
        auto inverseTrig = [&](const std::string& name, auto fn) {
            addFn(name, [name, fn](const std::vector<Value>& a)->Value {
                requireCount(a, 1, name); auto q = qarg(a, 0, name); requireDimensionless(q, name);
                return Quantity{ fn(q.value) };
                });
            };
        inverseTrig("asin", static_cast<double(*)(double)>(std::asin));
        inverseTrig("acos", static_cast<double(*)(double)>(std::acos));
        inverseTrig("atan", static_cast<double(*)(double)>(std::atan));
        addFn("atan2", [](const std::vector<Value>& a)->Value {
            requireCount(a, 2, "atan2"); auto x = qarg(a, 0, "atan2"), y = qarg(a, 1, "atan2");
            if (!(x.dim == y.dim)) throw ApsisError("atan2 arguments must have matching dimensions");
            return Quantity{ std::atan2(x.value,y.value) };
            });
        addFn("abs", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "abs");
            if (auto q = std::get_if<Quantity>(&a[0])) return Quantity{ std::abs(q->value),q->dim };
            if (auto v = std::get_if<Vec3>(&a[0])) return v->norm();
            throw ApsisError("abs requires a quantity or vector");
            });
        addFn("exp", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "exp"); auto q = qarg(a, 0, "exp"); requireDimensionless(q, "exp");
            return Quantity{ std::exp(q.value) };
            });
        addFn("log", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "log"); auto q = qarg(a, 0, "log"); requireDimensionless(q, "log");
            if (q.value <= 0) throw ApsisError("log requires a positive value");
            return Quantity{ std::log(q.value) };
            });
        addFn("log10", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "log10"); auto q = qarg(a, 0, "log10"); requireDimensionless(q, "log10");
            if (q.value <= 0) throw ApsisError("log10 requires a positive value");
            return Quantity{ std::log10(q.value) };
            });
        addFn("floor", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "floor"); auto q = qarg(a, 0, "floor"); return Quantity{ std::floor(q.value),q.dim };
            });
        addFn("ceil", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "ceil"); auto q = qarg(a, 0, "ceil"); return Quantity{ std::ceil(q.value),q.dim };
            });
        addFn("vec", [](const std::vector<Value>& a)->Value {
            requireCount(a, 3, "vec"); return Vec3{ qarg(a,0,"vec"),qarg(a,1,"vec"),qarg(a,2,"vec") };
            });
        addFn("norm", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "norm"); if (auto v = std::get_if<Vec3>(&a[0])) return v->norm();
            throw ApsisError("norm requires a vector");
            });
        addFn("unit", [](const std::vector<Value>& a)->Value {
            requireCount(a, 1, "unit"); if (auto v = std::get_if<Vec3>(&a[0])) return v->unit();
            throw ApsisError("unit requires a vector");
            });
        addFn("dot", [](const std::vector<Value>& a)->Value {
            requireCount(a, 2, "dot");
            auto x = std::get_if<Vec3>(&a[0]); auto y = std::get_if<Vec3>(&a[1]);
            if (!x || !y) throw ApsisError("dot requires two vectors");
            return x->dot(*y);
            });
        addFn("cross", [](const std::vector<Value>& a)->Value {
            requireCount(a, 2, "cross");
            auto x = std::get_if<Vec3>(&a[0]); auto y = std::get_if<Vec3>(&a[1]);
            if (!x || !y) throw ApsisError("cross requires two vectors");
            return x->cross(*y);
            });
        addFn("print", [](const std::vector<Value>& a)->Value {
            for (std::size_t i = 0;i < a.size();++i) { if (i) std::cout << ' '; std::cout << valueText(a[i]); }
            std::cout << '\n'; return {};
            });
    }

    static void installPackage(Environment& env, const std::string& packageName) {
        auto addFnIfMissing = [&](const std::string& name,
            std::function<Value(const std::vector<Value>&)> fn) {
                if (!env.values.contains(name))
                    env.define(name, nativeFunction(std::move(fn)), true);
            };
        if (packageName == "ephemerides") {
            addFnIfMissing("to_julian_date", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 6, "to_julian_date");
                auto integerArg = [&](std::size_t i) {
                    double v = qarg(a, i, "to_julian_date").value;
                    if (!std::isfinite(v) || std::floor(v) != v)
                        throw ApsisError("to_julian_date year, month, and day must be integers");
                    return static_cast<int>(v);
                    };
                return Quantity{ ApsisPackages::Ephemerides::to_julian_date(
                    integerArg(0), integerArg(1), integerArg(2),
                    qarg(a, 3, "to_julian_date").value,
                    qarg(a, 4, "to_julian_date").value,
                    qarg(a, 5, "to_julian_date").value) };
                });
            addFnIfMissing("centuries_from_j2000", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "centuries_from_j2000");
                return Quantity{ ApsisPackages::Ephemerides::centuries_from_j2000(qarg(a, 0, "centuries_from_j2000").value) };
                });
            addFnIfMissing("greenwich_mean_sidereal_time", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "greenwich_mean_sidereal_time");
                return Quantity{ ApsisPackages::Ephemerides::greenwich_mean_sidereal_time(qarg(a, 0, "greenwich_mean_sidereal_time").value) };
                });
            addFnIfMissing("get_sun_position_ecliptic", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "get_sun_position_ecliptic");
                double x, y, z;
                ApsisPackages::Ephemerides::get_sun_position_ecliptic(qarg(a, 0, "get_sun_position_ecliptic").value, x, y, z);
                Dimension length; length.e = { 1,0,0,0,0 };
                return Vec3{ {x,length},{y,length},{z,length} };
                });
            auto registerTransform = [&](const std::string& name, auto fn) {
                addFnIfMissing(name, [name, fn](const std::vector<Value>& a) -> Value {
                    requireCount(a, 4, name);
                    double x, y, z;
                    fn(qarg(a, 0, name).value, qarg(a, 1, name).value, qarg(a, 2, name).value,
                        qarg(a, 3, name).value, x, y, z);
                    Dimension length; length.e = { 1,0,0,0,0 };
                    return Vec3{ {x,length},{y,length},{z,length} };
                    });
                };
            registerTransform("ecliptic_to_equatorial", ApsisPackages::Ephemerides::ecliptic_to_equatorial);
            registerTransform("eci_to_ecef", ApsisPackages::Ephemerides::eci_to_ecef);
            registerTransform("ecef_to_eci", ApsisPackages::Ephemerides::ecef_to_eci);
            return;
        }

        if (packageName == "propagator") {
            addFnIfMissing("compute_density", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "compute_density");
                return Quantity{ ApsisPackages::Propagator::compute_density(qarg(a, 0, "compute_density").value) };
                });
            addFnIfMissing("propagate_rk4_step_position", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 8, "propagate_rk4_step_position");
                ApsisPackages::Propagator::State3D state{
                    qarg(a,0,"propagate_rk4_step_position").value,
                    qarg(a,1,"propagate_rk4_step_position").value,
                    qarg(a,2,"propagate_rk4_step_position").value,
                    qarg(a,3,"propagate_rk4_step_position").value,
                    qarg(a,4,"propagate_rk4_step_position").value,
                    qarg(a,5,"propagate_rk4_step_position").value
                };
                auto next = ApsisPackages::Propagator::propagate_rk4_step(
                    state, qarg(a, 6, "propagate_rk4_step_position").value,
                    qarg(a, 7, "propagate_rk4_step_position").value);
                Dimension length; length.e = { 1,0,0,0,0 };
                return Vec3{ {next.x,length},{next.y,length},{next.z,length} };
                });
            return;
        }

        if (packageName == "trajectory") {
            addFnIfMissing("stumpff_c2", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "stumpff_c2");
                return Quantity{ ApsisPackages::Trajectory::stumpff_c2(qarg(a, 0, "stumpff_c2").value) };
                });
            addFnIfMissing("stumpff_c3", [](const std::vector<Value>& a) -> Value {
                requireCount(a, 1, "stumpff_c3");
                return Quantity{ ApsisPackages::Trajectory::stumpff_c3(qarg(a, 0, "stumpff_c3").value) };
                });
            return;
        }

        if (packageName != "formulas") {
            throw ApsisError("unknown package '" + packageName + "' (available packages: formulas, ephemerides, propagator, trajectory)");
        }

        // Importing the same package more than once is harmless.
        if (env.values.contains("projectile_max_heigth")) return;

        auto addConstantIfMissing = [&](const std::string& name, double value, Dimension d = {}) {
            if (!env.values.contains(name)) env.define(name, Quantity{ value, d }, true);
            };
        auto dim = [](int l, int m, int t, int i = 0, int k = 0) {
            Dimension d; d.e = { l, m, t, i, k }; return d;
            };

        // Public constants exported by the formulas package. Existing Apsis constants
        // such as PI and G are kept, so importing does not replace built-ins.
        addConstantIfMissing("C", 2.99792458e8, dim(1, 0, -1));
        addConstantIfMissing("H", 6.62607015e-34, dim(2, 1, -1));
        addConstantIfMissing("H_BAR", 1.054571817e-34, dim(2, 1, -1));
        addConstantIfMissing("SIGMA", 5.670374e-8, dim(0, 1, -3, 0, -4));
        addConstantIfMissing("K_B", 1.380649e-23, dim(2, 1, -2, 0, -1));
        addConstantIfMissing("R_GAS", 8.31446261815324, dim(2, 1, -2, 0, -1));
        addConstantIfMissing("G0", 9.80665, dim(1, 0, -2));
        addConstantIfMissing("EPSILON_0", 8.8541878128e-12, dim(-3, -1, 4, 2));
        addConstantIfMissing("MU_0", 1.25663706e-6, dim(1, 1, -2, -2));
        addConstantIfMissing("K_E", 8.9875517923e9, dim(3, 1, -4, -2));
        addConstantIfMissing("WIEN_B", 2.897771955e-3, dim(1, 0, 0, 0, 1));
        addConstantIfMissing("M_P", 1.67262192369e-27, dim(0, 1, 0));
        addConstantIfMissing("M_H", 1.6735575e-27, dim(0, 1, 0));
        addConstantIfMissing("SIGMA_T", 6.6524587321e-29, dim(2, 0, 0));

        // Formula contracts: dimensions are (length, mass, time, current, temperature).
        // Bare numeric inputs remain accepted as SI-valued numbers for compatibility;
        // explicitly dimensioned inputs must match the declared dimensions.
        const Dimension Z = dim(0, 0, 0), L = dim(1, 0, 0), M = dim(0, 1, 0),
            T = dim(0, 0, 1), I = dim(0, 0, 0, 1), K = dim(0, 0, 0, 0, 1),
            V = dim(1, 0, -1), A = dim(1, 0, -2), F = dim(1, 1, -2),
            E = dim(2, 1, -2), P = dim(2, 1, -3), PR = dim(-1, 1, -2),
            DENS = dim(-3, 1, 0), MOM = dim(1, 1, -1), IMP = dim(1, 1, -1),
            MOI = dim(2, 1, 0), AV = dim(0, 0, -1), AA = dim(0, 0, -2),
            CHARGE = dim(0, 0, 1, 1), VOLT = dim(2, 1, -3, -1),
            RESIST = dim(2, 1, -3, -2), CAP = dim(-2, -1, 4, 2),
            INDUCT = dim(2, 1, -2, -2), FLUX = dim(2, 1, -2, -1),
            EFIELD = dim(1, 1, -3, -1), BFIELD = dim(0, 1, -2, -1),
            INTENSITY = dim(0, 1, -3), ENTROPY = dim(2, 1, -2, 0, -1),
            SPECHEAT = dim(2, 0, -2, 0, -1), CONDUCT = dim(1, 1, -3, 0, -1),
            DYNAMICVISC = dim(-1, 1, -1), STRESS = PR, STRAIN = Z;

        struct FormulaSignature { std::vector<Dimension> inputs; Dimension output; };
        std::unordered_map<std::string, FormulaSignature> signatures;
        auto addSig = [&](const std::string& n, std::initializer_list<Dimension> in, Dimension out) {
            signatures.emplace(n, FormulaSignature{ std::vector<Dimension>(in), out });
            };

        addSig("projectile_max_heigth", { V,Z,A }, L);
        addSig("projectile_range", { V,Z,A }, L);
        addSig("projectile_time_of_flight", { V,Z,A }, T);
        addSig("average_velocity", { L,T }, V);
        addSig("average_acceleration", { V,T }, A);
        addSig("kinematic_final_velocity", { V,A,T }, V);
        addSig("kinematic_position", { L,V,A,T }, L);
        addSig("kinematic_velocity_squared", { V,A,L }, dim(2, 0, -2));
        addSig("newton_second_law", { M,A }, F);
        addSig("weight", { M,A }, F);
        addSig("gravitational_force", { M,M,L }, F);
        addSig("linear_momentum", { M,V }, MOM);
        addSig("impulse", { F,T }, IMP);
        addSig("friction_force", { Z,F }, F);
        addSig("work_done", { F,L,Z }, E);
        addSig("kinetic_energy", { M,V }, E);
        addSig("gravitational_potential_energy", { M,A,L }, E);
        addSig("mechanical_power", { E,T }, P);
        addSig("hookes_law", { dim(0,1,-2),L }, F);
        addSig("elastic_potential_energy", { dim(0,1,-2),L }, E);
        addSig("angular_velocity", { Z,T }, AV);
        addSig("angular_acceleration", { AV,T }, AA);
        addSig("tangential_velocity", { AV,L }, V);
        addSig("centripetal_acceleration", { V,L }, A);
        addSig("centripetal_force", { M,V,L }, F);
        addSig("torque", { L,F,Z }, E);
        addSig("moment_of_inertia_point_mass", { M,L }, MOI);
        addSig("moment_of_inertia_solid_disk", { M,L }, MOI);
        addSig("moment_of_inertia_solid_sphere", { M,L }, MOI);
        addSig("moment_of_inertia_thin_rod_center", { M,L }, MOI);
        addSig("angular_momentum", { MOI,AV }, dim(2, 1, -1));
        addSig("rotational_kinetic_energy", { MOI,AV }, E);
        addSig("parallel_axis_theorem", { MOI,M,L }, MOI);
        addSig("shm_spring_period", { M,dim(0,1,-2) }, T);
        addSig("simple_pendulum_period", { L,A }, T);
        addSig("shm_position", { L,AV,T,Z }, L);
        addSig("shm_angular_frequency", { dim(0,1,-2),M }, AV);
        addSig("damped_oscillation_amplitude", { L,dim(0,1,-1),M,T }, L);
        addSig("ideal_gas_pressure", { Z,K,L.scaled(3) }, PR);
        addSig("first_law_thermodynamics", { E,E }, E);
        addSig("sensible_heat", { M,SPECHEAT,K }, E);
        addSig("latent_heat", { M,dim(2,0,-2) }, E);
        addSig("carnot_efficiency", { K,K }, Z);
        addSig("entropy_change_reversible", { E,K }, ENTROPY);
        addSig("statistical_entropy", { Z }, dim(2, 1, -2, 0, -1));
        addSig("stefan_boltzmann_power", { L.scaled(2),K }, P);
        addSig("wiens_displacement_law", { K }, L);
        addSig("linear_thermal_expansion", { dim(0,0,0,0,-1),L,K }, L);
        addSig("fourier_conduction_rate", { CONDUCT,L.scaled(2),K,L }, P);
        addSig("rms_molecular_speed", { K,dim(0,1,0) }, V);
        addSig("isobaric_work", { PR,L.scaled(3) }, E);
        addSig("coulomb_force", { CHARGE,CHARGE,L }, F);
        addSig("electric_field_point_charge", { CHARGE,L }, EFIELD);
        addSig("electric_force_on_charge", { CHARGE,EFIELD }, F);
        addSig("electric_potential_point_charge", { CHARGE,L }, VOLT);
        addSig("parallel_plate_capacitance", { L.scaled(2),L }, CAP);
        addSig("capacitor_stored_energy", { CAP,VOLT }, E);
        addSig("ohms_law_voltage", { I,RESIST }, VOLT);
        addSig("electrical_power", { VOLT,I }, P);
        addSig("conductor_resistance", { dim(3,1,-3,-2),L,L.scaled(2) }, RESIST);
        addSig("amperes_law_magnetic_field_circulation", { I }, BFIELD);
        addSig("faradays_law_emf", { FLUX,T }, VOLT);
        addSig("self_inductance_emf", { INDUCT,I,T }, VOLT);
        addSig("inductor_stored_energy", { INDUCT,I }, E);
        addSig("speed_of_em_waves", {}, V);
        addSig("wave_speed", { dim(0,0,-1),L }, V);
        addSig("wave_period", { dim(0,0,-1) }, T);
        addSig("snells_law_angle2", { Z,Z,Z }, Z);
        addSig("critical_angle", { Z,Z }, Z);
        addSig("doppler_effect_observed_frequency", { dim(0,0,-1),V,V,V }, dim(0, 0, -1));
        addSig("thin_lens_image_distance", { L,L }, L);
        addSig("lens_magnification", { L,L }, Z);
        addSig("wave_intensity", { P,L.scaled(2) }, INTENSITY);
        addSig("sound_level_decibels", { INTENSITY }, Z);
        addSig("mass_energy_equivalence", { M }, E);
        addSig("lorentz_factor", { V }, Z);
        addSig("time_dilation", { T,Z }, T);
        addSig("length_contraction", { L,Z }, L);
        addSig("relativistic_momentum", { M,V,Z }, MOM);
        addSig("relativistic_total_energy", { M,Z }, E);
        addSig("energy_momentum_relation", { MOM,M }, E);
        addSig("photon_energy_frequency", { dim(0,0,-1) }, E);
        addSig("de_broglie_wavelength", { MOM }, L);
        addSig("radioactive_decay_remaining", { Z,dim(0,0,-1),T }, Z);
        addSig("radioactive_half_life", { dim(0,0,-1) }, T);
        addSig("keplers_third_law", { L,M }, T);
        addSig("stellar_luminosity", { L,K }, P);
        addSig("flux_luminosity_relation", { P,L }, INTENSITY);
        addSig("distance_modulus", { L }, Z);
        addSig("parallax_distance", { Z }, L);
        addSig("hubbles_law_velocity", { dim(0,0,-1),L }, V);
        addSig("cosmological_redshift", { L,L }, Z);
        addSig("schwarzschild_radius", { M }, L);
        addSig("eddington_luminosity", { M }, P);
        addSig("vis_viva_velocity", { dim(3,-1,-2),M,L,L }, V);
        addSig("circular_orbital_velocity", { M,L }, V);
        addSig("escape_velocity", { M,L }, V);
        addSig("specific_orbital_energy", { V,M,L }, dim(2, 0, -2));
        addSig("orbital_eccentricity_from_apsides", { L,L }, Z);
        addSig("semi_major_axis_from_apsides", { L,L }, L);
        addSig("hohmann_transfer_delta_v1", { M,L,L }, V);
        addSig("hohmann_transfer_delta_v2", { M,L,L }, V);
        addSig("sphere_of_influence_radius", { L,M,M }, L);
        addSig("tsiolkovsky_rocket_equation", { V,M,M }, V);
        addSig("effective_exhaust_velocity", { T,A }, V);
        addSig("specific_impulse", { F,dim(0,1,-1),A }, T);
        addSig("thrust_equation", { dim(0,1,-1),V,PR,PR,L.scaled(2) }, F);
        addSig("rocket_mass_ratio", { M,M }, Z);
        addSig("propellant_mass_fraction", { M,M }, Z);
        addSig("thrust_to_weight_ratio", { F,M,A }, Z);
        addSig("rocket_burn_time", { M,dim(0,1,-1) }, T);
        addSig("normal_stress", { F,L.scaled(2) }, STRESS);
        addSig("normal_strain", { L,L }, STRAIN);
        addSig("youngs_modulus", { STRESS,STRAIN }, PR);
        addSig("shear_modulus", { STRESS,STRAIN }, PR);
        addSig("poissons_ratio", { Z,Z }, Z);
        addSig("bending_stress", { E,L,MOI }, STRESS);
        addSig("torsional_shear_stress", { E,L,L.scaled(4) }, STRESS);
        addSig("angle_of_twist", { E,L,L.scaled(4),PR }, Z);
        addSig("euler_buckling_load", { PR,L.scaled(4),Z,L }, F);
        addSig("factor_of_safety", { STRESS,STRESS }, Z);
        addSig("thermal_stress_constrained", { PR,dim(0,0,0,0,-1),K }, STRESS);
        addSig("cantilever_beam_deflection", { F,L,PR,L.scaled(4) }, L);
        addSig("continuity_equation_v2", { L.scaled(2),V,L.scaled(2) }, V);
        addSig("reynolds_number", { DENS,V,L,DYNAMICVISC }, Z);
        addSig("hydrostatic_pressure", { PR,DENS,A,L }, PR);
        addSig("buoyant_force", { DENS,A,L.scaled(3) }, F);
        addSig("drag_force", { Z,DENS,L.scaled(2),V }, F);
        addSig("lift_force", { Z,DENS,L.scaled(2),V }, F);
        addSig("gear_ratio", { Z,Z }, Z);
        addSig("newtons_law_cooling_convection", { dim(1,1,-3,0,-1),L.scaled(2),K }, P);

        auto addFormula = [&](const std::string& name, std::size_t count,
            std::function<double(const std::vector<double>&)> fn) {
                auto sigIt = signatures.find(name);
                if (sigIt == signatures.end())
                    throw ApsisError("missing dimensional signature for formula '" + name + "'");
                const FormulaSignature signature = sigIt->second;
                if (signature.inputs.size() != count)
                    throw ApsisError("internal signature arity mismatch for formula '" + name + "'");
                env.define(name, nativeFunction([name, count, signature, fn = std::move(fn)](const std::vector<Value>& args) -> Value {
                    requireCount(args, count, name);
                    std::vector<double> numbers;
                    numbers.reserve(args.size());
                    for (std::size_t i = 0; i < args.size(); ++i) {
                        Quantity q = asQuantity(args[i], name + " argument " + std::to_string(i + 1));
                        const Dimension expected = signature.inputs[i];
                        // Bare numbers are SI values; unit-bearing values are checked strictly.
                        if (!q.dim.dimensionless() && !(q.dim == expected)) {
                            throw ApsisError(name + " argument " + std::to_string(i + 1) +
                                " has dimension " + dimensionName(q.dim) +
                                "; expected " + dimensionName(expected));
                        }
                        // The legacy astronomy helpers use parsecs/arcseconds as bare numbers.
                        // Explicit length quantities are SI metres, so convert metres to parsecs.
                        constexpr double PARSEC_METRES = 3.0856775814913673e16;
                        if (name == "distance_modulus" && i == 0 && !q.dim.dimensionless())
                            numbers.push_back(q.value / PARSEC_METRES);
                        else
                            numbers.push_back(q.value);
                    }
                    try {
                        double result = fn(numbers);
                        // parallax_distance returns a legacy numeric parsec count; the runtime
                        // represents physical lengths in SI metres.
                        if (name == "parallax_distance") result *= 3.0856775814913673e16;
                        return Quantity{ result, signature.output };
                    }
                    catch (const ApsisError&) { throw; }
                    catch (const std::exception& ex) { throw ApsisError(ex.what()); }
                    }), true);
            };

        addFormula("projectile_max_heigth", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::projectile_max_heigth(x[0], x[1], x[2]); });
        addFormula("projectile_range", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::projectile_range(x[0], x[1], x[2]); });
        addFormula("projectile_time_of_flight", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::projectile_time_of_flight(x[0], x[1], x[2]); });
        addFormula("average_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::average_velocity(x[0], x[1]); });
        addFormula("average_acceleration", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::average_acceleration(x[0], x[1]); });
        addFormula("kinematic_final_velocity", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::kinematic_final_velocity(x[0], x[1], x[2]); });
        addFormula("kinematic_position", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::kinematic_position(x[0], x[1], x[2], x[3]); });
        addFormula("kinematic_velocity_squared", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::kinematic_velocity_squared(x[0], x[1], x[2]); });
        addFormula("newton_second_law", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::newton_second_law(x[0], x[1]); });
        addFormula("weight", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::weight(x[0], x[1]); });
        addFormula("gravitational_force", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::gravitational_force(x[0], x[1], x[2]); });
        addFormula("linear_momentum", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::linear_momentum(x[0], x[1]); });
        addFormula("impulse", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::impulse(x[0], x[1]); });
        addFormula("friction_force", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::friction_force(x[0], x[1]); });
        addFormula("work_done", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::work_done(x[0], x[1], x[2]); });
        addFormula("kinetic_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::kinetic_energy(x[0], x[1]); });
        addFormula("gravitational_potential_energy", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::gravitational_potential_energy(x[0], x[1], x[2]); });
        addFormula("mechanical_power", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::mechanical_power(x[0], x[1]); });
        addFormula("hookes_law", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::hookes_law(x[0], x[1]); });
        addFormula("elastic_potential_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::elastic_potential_energy(x[0], x[1]); });
        addFormula("angular_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::angular_velocity(x[0], x[1]); });
        addFormula("angular_acceleration", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::angular_acceleration(x[0], x[1]); });
        addFormula("tangential_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::tangential_velocity(x[0], x[1]); });
        addFormula("centripetal_acceleration", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::centripetal_acceleration(x[0], x[1]); });
        addFormula("centripetal_force", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::centripetal_force(x[0], x[1], x[2]); });
        addFormula("torque", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::torque(x[0], x[1], x[2]); });
        addFormula("moment_of_inertia_point_mass", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::moment_of_inertia_point_mass(x[0], x[1]); });
        addFormula("moment_of_inertia_solid_disk", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::moment_of_inertia_solid_disk(x[0], x[1]); });
        addFormula("moment_of_inertia_solid_sphere", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::moment_of_inertia_solid_sphere(x[0], x[1]); });
        addFormula("moment_of_inertia_thin_rod_center", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::moment_of_inertia_thin_rod_center(x[0], x[1]); });
        addFormula("angular_momentum", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::angular_momentum(x[0], x[1]); });
        addFormula("rotational_kinetic_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::rotational_kinetic_energy(x[0], x[1]); });
        addFormula("parallel_axis_theorem", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::parallel_axis_theorem(x[0], x[1], x[2]); });
        addFormula("shm_spring_period", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::shm_spring_period(x[0], x[1]); });
        addFormula("simple_pendulum_period", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::simple_pendulum_period(x[0], x[1]); });
        addFormula("shm_position", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::shm_position(x[0], x[1], x[2], x[3]); });
        addFormula("shm_angular_frequency", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::shm_angular_frequency(x[0], x[1]); });
        addFormula("damped_oscillation_amplitude", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::damped_oscillation_amplitude(x[0], x[1], x[2], x[3]); });
        addFormula("ideal_gas_pressure", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::ideal_gas_pressure(x[0], x[1], x[2]); });
        addFormula("first_law_thermodynamics", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::first_law_thermodynamics(x[0], x[1]); });
        addFormula("sensible_heat", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::sensible_heat(x[0], x[1], x[2]); });
        addFormula("latent_heat", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::latent_heat(x[0], x[1]); });
        addFormula("carnot_efficiency", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::carnot_efficiency(x[0], x[1]); });
        addFormula("entropy_change_reversible", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::entropy_change_reversible(x[0], x[1]); });
        addFormula("statistical_entropy", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::statistical_entropy(x[0]); });
        addFormula("stefan_boltzmann_power", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::stefan_boltzmann_power(x[0], x[1]); });
        addFormula("wiens_displacement_law", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::wiens_displacement_law(x[0]); });
        addFormula("linear_thermal_expansion", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::linear_thermal_expansion(x[0], x[1], x[2]); });
        addFormula("fourier_conduction_rate", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::fourier_conduction_rate(x[0], x[1], x[2], x[3]); });
        addFormula("rms_molecular_speed", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::rms_molecular_speed(x[0], x[1]); });
        addFormula("isobaric_work", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::isobaric_work(x[0], x[1]); });
        addFormula("coulomb_force", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::coulomb_force(x[0], x[1], x[2]); });
        addFormula("electric_field_point_charge", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::electric_field_point_charge(x[0], x[1]); });
        addFormula("electric_force_on_charge", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::electric_force_on_charge(x[0], x[1]); });
        addFormula("electric_potential_point_charge", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::electric_potential_point_charge(x[0], x[1]); });
        addFormula("parallel_plate_capacitance", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::parallel_plate_capacitance(x[0], x[1]); });
        addFormula("capacitor_stored_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::capacitor_stored_energy(x[0], x[1]); });
        addFormula("ohms_law_voltage", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::ohms_law_voltage(x[0], x[1]); });
        addFormula("electrical_power", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::electrical_power(x[0], x[1]); });
        addFormula("conductor_resistance", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::conductor_resistance(x[0], x[1], x[2]); });
        addFormula("amperes_law_magnetic_field_circulation", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::amperes_law_magnetic_field_circulation(x[0]); });
        addFormula("faradays_law_emf", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::faradays_law_emf(x[0], x[1]); });
        addFormula("self_inductance_emf", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::self_inductance_emf(x[0], x[1], x[2]); });
        addFormula("inductor_stored_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::inductor_stored_energy(x[0], x[1]); });
        addFormula("speed_of_em_waves", 0, [](const std::vector<double>&) { return ApsisPackages::Formulas::speed_of_em_waves(); });
        addFormula("wave_speed", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::wave_speed(x[0], x[1]); });
        addFormula("wave_period", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::wave_period(x[0]); });
        addFormula("snells_law_angle2", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::snells_law_angle2(x[0], x[1], x[2]); });
        addFormula("critical_angle", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::critical_angle(x[0], x[1]); });
        addFormula("doppler_effect_observed_frequency", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::doppler_effect_observed_frequency(x[0], x[1], x[2], x[3]); });
        addFormula("thin_lens_image_distance", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::thin_lens_image_distance(x[0], x[1]); });
        addFormula("lens_magnification", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::lens_magnification(x[0], x[1]); });
        addFormula("wave_intensity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::wave_intensity(x[0], x[1]); });
        addFormula("sound_level_decibels", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::sound_level_decibels(x[0]); });
        addFormula("mass_energy_equivalence", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::mass_energy_equivalence(x[0]); });
        addFormula("lorentz_factor", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::lorentz_factor(x[0]); });
        addFormula("time_dilation", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::time_dilation(x[0], x[1]); });
        addFormula("length_contraction", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::length_contraction(x[0], x[1]); });
        addFormula("relativistic_momentum", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::relativistic_momentum(x[0], x[1], x[2]); });
        addFormula("relativistic_total_energy", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::relativistic_total_energy(x[0], x[1]); });
        addFormula("energy_momentum_relation", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::energy_momentum_relation(x[0], x[1]); });
        addFormula("photon_energy_frequency", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::photon_energy_frequency(x[0]); });
        addFormula("de_broglie_wavelength", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::de_broglie_wavelength(x[0]); });
        addFormula("radioactive_decay_remaining", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::radioactive_decay_remaining(x[0], x[1], x[2]); });
        addFormula("radioactive_half_life", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::radioactive_half_life(x[0]); });
        addFormula("keplers_third_law", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::keplers_third_law(x[0], x[1]); });
        addFormula("stellar_luminosity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::stellar_luminosity(x[0], x[1]); });
        addFormula("flux_luminosity_relation", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::flux_luminosity_relation(x[0], x[1]); });
        addFormula("distance_modulus", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::distance_modulus(x[0]); });
        addFormula("parallax_distance", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::parallax_distance(x[0]); });
        addFormula("hubbles_law_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::hubbles_law_velocity(x[0], x[1]); });
        addFormula("cosmological_redshift", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::cosmological_redshift(x[0], x[1]); });
        addFormula("schwarzschild_radius", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::schwarzschild_radius(x[0]); });
        addFormula("eddington_luminosity", 1, [](const std::vector<double>& x) { return ApsisPackages::Formulas::eddington_luminosity(x[0]); });
        addFormula("vis_viva_velocity", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::vis_viva_velocity(x[0], x[1], x[2], x[3]); });
        addFormula("circular_orbital_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::circular_orbital_velocity(x[0], x[1]); });
        addFormula("escape_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::escape_velocity(x[0], x[1]); });
        addFormula("specific_orbital_energy", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::specific_orbital_energy(x[0], x[1], x[2]); });
        addFormula("orbital_eccentricity_from_apsides", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::orbital_eccentricity_from_apsides(x[0], x[1]); });
        addFormula("semi_major_axis_from_apsides", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::semi_major_axis_from_apsides(x[0], x[1]); });
        addFormula("hohmann_transfer_delta_v1", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::hohmann_transfer_delta_v1(x[0], x[1], x[2]); });
        addFormula("hohmann_transfer_delta_v2", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::hohmann_transfer_delta_v2(x[0], x[1], x[2]); });
        addFormula("sphere_of_influence_radius", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::sphere_of_influence_radius(x[0], x[1], x[2]); });
        addFormula("tsiolkovsky_rocket_equation", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::tsiolkovsky_rocket_equation(x[0], x[1], x[2]); });
        addFormula("effective_exhaust_velocity", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::effective_exhaust_velocity(x[0], x[1]); });
        addFormula("specific_impulse", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::specific_impulse(x[0], x[1], x[2]); });
        addFormula("thrust_equation", 5, [](const std::vector<double>& x) { return ApsisPackages::Formulas::thrust_equation(x[0], x[1], x[2], x[3], x[4]); });
        addFormula("rocket_mass_ratio", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::rocket_mass_ratio(x[0], x[1]); });
        addFormula("propellant_mass_fraction", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::propellant_mass_fraction(x[0], x[1]); });
        addFormula("thrust_to_weight_ratio", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::thrust_to_weight_ratio(x[0], x[1], x[2]); });
        addFormula("rocket_burn_time", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::rocket_burn_time(x[0], x[1]); });
        addFormula("normal_stress", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::normal_stress(x[0], x[1]); });
        addFormula("normal_strain", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::normal_strain(x[0], x[1]); });
        addFormula("youngs_modulus", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::youngs_modulus(x[0], x[1]); });
        addFormula("shear_modulus", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::shear_modulus(x[0], x[1]); });
        addFormula("poissons_ratio", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::poissons_ratio(x[0], x[1]); });
        addFormula("bending_stress", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::bending_stress(x[0], x[1], x[2]); });
        addFormula("torsional_shear_stress", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::torsional_shear_stress(x[0], x[1], x[2]); });
        addFormula("angle_of_twist", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::angle_of_twist(x[0], x[1], x[2], x[3]); });
        addFormula("euler_buckling_load", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::euler_buckling_load(x[0], x[1], x[2], x[3]); });
        addFormula("factor_of_safety", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::factor_of_safety(x[0], x[1]); });
        addFormula("thermal_stress_constrained", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::thermal_stress_constrained(x[0], x[1], x[2]); });
        addFormula("cantilever_beam_deflection", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::cantilever_beam_deflection(x[0], x[1], x[2], x[3]); });
        addFormula("continuity_equation_v2", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::continuity_equation_v2(x[0], x[1], x[2]); });
        addFormula("reynolds_number", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::reynolds_number(x[0], x[1], x[2], x[3]); });
        addFormula("hydrostatic_pressure", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::hydrostatic_pressure(x[0], x[1], x[2], x[3]); });
        addFormula("buoyant_force", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::buoyant_force(x[0], x[1], x[2]); });
        addFormula("drag_force", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::drag_force(x[0], x[1], x[2], x[3]); });
        addFormula("lift_force", 4, [](const std::vector<double>& x) { return ApsisPackages::Formulas::lift_force(x[0], x[1], x[2], x[3]); });
        addFormula("gear_ratio", 2, [](const std::vector<double>& x) { return ApsisPackages::Formulas::gear_ratio(x[0], x[1]); });
        addFormula("newtons_law_cooling_convection", 3, [](const std::vector<double>& x) { return ApsisPackages::Formulas::newtons_law_cooling_convection(x[0], x[1], x[2]); });
    }

    static std::string readFile(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) throw ApsisError("could not open script: " + path);
        std::ostringstream contents;
        contents << in.rdbuf();
        return contents.str();
    }

    static int run(const std::string& source, const std::string& label) {
        try {
            Lexer lexer(source);
            Parser parser(lexer.scan());
            auto statements = parser.parse();
            auto global = std::make_shared<Environment>();
            installBuiltins(global);
            for (auto& s : statements) s->exec(*global);
            return 0;
        }
        catch (const ReturnSignal&) {
            std::cerr << label << ": return used outside a function\n";
            return 1;
        }
        catch (const ApsisError& ex) {
            std::cerr << label << ": " << ex.what() << '\n';
            return 1;
        }
        catch (const std::exception& ex) {
            std::cerr << label << ": internal error: " << ex.what() << '\n';
            return 1;
        }
    }

    int runFile(const std::string& path) {
        try {
            return run(readFile(path), path);
        }
        catch (const ApsisError& ex) {
            std::cerr << "Apsis: " << ex.what() << '\n';
            return 1;
        }
    }

}
