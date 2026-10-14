# Apsis implementation notes

## Script language
The `.jdx` interpreter retains the existing `const`, `define`, `fn`, `return`, `if`/`else`, `while`, `print`, `import formulas`, quantity arithmetic, vectors and dimensional signatures.

## Numerical API
`apsis::numerics` provides bracketed bisection, Newton iteration with central finite-difference derivative, adaptive Simpson quadrature, fixed-step RK4, and adaptive Dormand–Prince 5(4). Numerical methods reject invalid tolerances, non-finite results and exhausted iteration/step budgets. The adaptive RK45 API returns the final state and accepted/rejected step counts; it does not yet provide dense output or event root isolation.

## Orbital API
`apsis::orbit` uses SI units: metres, metres/second, seconds and m^3/s^2. `Body::mu` is the central body's gravitational parameter. `propagateRK4` is a fixed-step Cartesian propagator. Optional force models include J2, third-body differential gravity, simplified atmospheric drag, constant-direction solar radiation pressure and thrust. `elementsFromState` returns classical osculating elements; `stateFromElements` currently supports elliptic orbits (`0 <= e < 1`).

## Formula package conventions
The legacy formula package accepts plain doubles. Most functions expect SI values, but astronomy helpers retain historical conventions: `distance_modulus` expects distance in parsecs and `parallax_distance` accepts parallax in arcseconds; the interpreter converts its parsec result to SI metres. For backward compatibility, bare numeric input to `distance_modulus` is interpreted as parsecs, while an explicit length quantity is converted from metres to parsecs. Angle arguments to formula functions are radians unless otherwise stated.

## Verification policy
A passing test verifies only the tested behavior and tolerance. Formula registration and dimensional signatures do not by themselves prove physical correctness. Reference comparisons, provenance, assumptions and domain limits must be documented for each scientific model before high-fidelity use.
