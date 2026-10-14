#pragma once
#include <functional>
#include <vector>
#include <cstddef>
#include <limits>

namespace apsis::numerics {
struct RootResult { double root; std::size_t iterations; double residual; };
RootResult bisection(const std::function<double(double)>& f, double a, double b,
                     double abs_tol=1e-12, std::size_t max_iter=200);
RootResult newton(const std::function<double(double)>& f, double x0,
                  double abs_tol=1e-12, std::size_t max_iter=100,
                  double derivative_step=1e-6);
double derivative(const std::function<double(double)>& f, double x, double h=1e-5);
double integrateAdaptiveSimpson(const std::function<double(double)>& f, double a, double b,
                                double abs_tol=1e-10, std::size_t max_depth=20);
using State = std::vector<double>;
using OdeFunction = std::function<State(double, const State&)>;
struct OdeResult { double t; State y; std::size_t accepted_steps; std::size_t rejected_steps; };
OdeResult integrateRK4(const OdeFunction& f, double t0, const State& y0, double t1, double step);
OdeResult integrateRK45(const OdeFunction& f, double t0, const State& y0, double t1,
                        double initial_step=0.1, double abs_tol=1e-9,
                        double rel_tol=1e-9, double min_step=1e-12,
                        double max_step=std::numeric_limits<double>::infinity());
}
