# stateEstimator
This repository contains a C implementation of an Extended Kalman Filter (EKF) and supporting matrix math utilities,
plus a suite of tests and tooling for documentation, formatting, and coverage. The core library lives in `src/` and
`inc/`, with focused unit/integration tests under `test/`. It is intended for embedded/robotics-style sensor fusion
workflows and includes examples that exercise linear and nonlinear models.

Key folders:
- `inc/`: Public headers for the EKF and matrix utilities.
- `src/`: Implementation sources.
- `test/`: Test binaries for EKF behavior and matrix math.
- `docs/`: Generated Doxygen output.

## Notes
- `EKFPredict` propagates state using `f(x)` and uses `updateAMatrix`/`A` for linear models or to support Jacobian generation.
- `EKFUpdate` uses the Joseph form covariance update for better numerical stability.
- If the innovation matrix `S` is near-singular, the update applies diagonal jitter to attempt inversion.
- `STATIC_MATRIX_DIRECTIVE` uses static storage, so it is not thread-safe and should not be shared across concurrent contexts.

## Documentation
To generate API docs with Doxygen:
1. Install `doxygen`.
2. From the repo root, run:
```powershell
doxygen Doxyfile
```
The HTML output will be under `docs/html`.

## Code Coverage
Coverage is supported for GCC/Clang builds using `gcovr`.
1. Configure with coverage enabled:
```powershell
cmake -S . -B build -DENABLE_COVERAGE=ON
```
2. Build and run the coverage target:
```powershell
cmake --build build
cmake --build build --target coverage
```
The reports will be written to `build/coverage.html` and `build/coverage.xml`.
