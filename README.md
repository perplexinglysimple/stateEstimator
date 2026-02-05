# stateEstimator
Takes an IMU, MPU-9250, and a XA1110 GPS and throws it in a EKF. This project will be written in C and run eventually on an arduino

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
