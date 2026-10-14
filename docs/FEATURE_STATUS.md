# Feature status (22-item checklist)

This package is a working modular baseline, not a claim that every roadmap item is finished.

| # | Area | Current implementation |
|---:|---|---|
| 1 | Automated tests and scientific validation | **Partial:** CTest covers numerical primitives, matrix solve, unit conversion, and orbital smoke cases. Not every formula has an independent reference test yet. |
| 2 | Complete SI units | **Partial:** standalone seven-base SI dimensions and conversion table; legacy interpreter still uses five dimensions. |
| 3 | Numerical computation | **Partial:** bisection, Newton, derivative, adaptive Simpson, RK4, RK45, dense matrices and linear solve. Brent, eigenvalues, least squares, and interpolation are not included yet. |
| 4 | ODE integration | **Partial:** fixed-step RK4 and adaptive Dormand–Prince 5(4) API. Dense output and event-coupled adaptive orbit integration remain missing. |
| 5 | Frames and reference systems | **Partial:** Cartesian state and frame labels/Z rotation; complete inertial/body-fixed velocity transformations and time scales remain missing. |
| 6 | Orbital propagation and perturbations | **Partial:** two-body, J2, third-body differential gravity, simplified drag/SRP, thrust, and classical elements. Higher gravity harmonics, finite-burn schedules and validated perturbation models remain missing. |
| 7 | Ephemerides and constants | **Partial:** explicit body `mu`, radius, J2, rotation-rate fields and CSV state import. No SPICE kernel/time conversion integration yet. |
| 8 | Events and mission analysis | **Partial:** analytical transfer/rocket formula package remains; robust trajectory event detection, close-approach screening and launch-window optimization remain missing. |
| 9 | Formula dimensions | **Partial:** formula signature checking is preserved; legacy bare numbers remain allowed for compatibility. Parsecs are converted to SI on the interpreter output path. Strict mode is not yet configurable. |
| 10 | Vector math | **Partial:** finite-safe vector normalization, dot/cross/norm and orbital vector operations. DSL vector API is not expanded in this revision. |
| 11 | Operator semantics | Existing parser preserved. A full specification-backed precedence/short-circuit regression suite is still needed. |
| 12 | Scope and control flow | Existing language semantics preserved. `break`/`continue` and complete stack diagnostics are not added. |
| 13 | Runtime safety | Numerical modules check non-finite values and impose step budgets. Configurable wall-clock, memory, source-size and call-stack budgets remain missing. |
| 14 | Formula registry | **Partial:** formula arity/signature registry is validated on import. References, provenance, domain metadata and per-formula validation cases remain missing. |
| 15 | Arrays | Not yet exposed in the `.jdx` language. |
| 16 | Matrices | **Partial:** standalone dense matrix multiply, transpose and Gaussian-elimination solve; not exposed in `.jdx`. |
| 17 | File I/O | **Partial:** state CSV import/export and XY CSV export; no general structured-data parser. |
| 18 | Modules/packages | **Partial:** existing `import formulas` preserved; no general package loader yet. |
| 19 | Plotting/export | **Partial:** CSV export for plotting in external tools; no native plotting UI. |
| 20 | Debugging/profiling | Not yet implemented beyond source-located parser/runtime diagnostics already present. |
| 21 | Specification/docs | **Partial:** README and science/language notes added; full formal grammar and semantics specification remain future work. |
| 22 | CI/portability | **Partial:** CMake and CTest work in this Linux build. Windows/macOS CI workflows are not configured. |
