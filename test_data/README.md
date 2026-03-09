# EKF Reference Validation

Python-based validation framework for comparing a target EKF output against a Python reference implementation.

## Purpose

This tool validates EKF behavior by:

1. Loading real sensor data from CSV files
2. Running a reference EKF implementation from a trusted python library with the same measurement sequence
3. Comparing state outputs and computing error metrics
4. Generating detailed reports with error analysis


The data was gathered from sitting on my bench. Assume a stationary system with no motion. The EKF should ideally output a constant state estimate with zero velocity and no change in orientation. Any deviation from this can be attributed to sensor noise, bias, or implementation errors.

We should model this as a "person". Allow the system model to walk around some. Future data will include more dynamic motion, but for now we want to validate the EKF's ability to handle a simple stationary scenario.

## Quick Start

This repo now includes a concrete C-vs-Python comparison path for the bench dataset:

1. Build the C export test.
2. Run the export test on the CSVs in this folder.
3. Replay the same data through a FilterPy `ExtendedKalmanFilter` reference.
4. Compare state and covariance traces row-by-row.

```bash
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
mkdir -p build
cd build
cmake ..
cmake --build . --target ekfBenchDataExportTest
./ekfBenchDataExportTest ekf_bench_compare.csv
../.venv/bin/python ../scripts/compare_filterpy_bench.py --csv ekf_bench_compare.csv
```

## Data Sources

- `ICM20948.CSV`: IMU accelerometer/gyro/magnetometer samples. The current comparison uses accelerometer data for prediction.
- `SAM_M8Q.CSV`: GPS fixes. The current comparison uses latitude/longitude when `valid_fix=1` and the coordinates are non-zero.
- `BMP384.CSV`: Barometer samples. The current comparison converts pressure to relative altitude.

## Current Model

The comparison uses a simple local-frame "person" model that is intentionally modest:

- State: `[x, y, z, vx, vy, vz]`
- Prediction: constant-velocity kinematics driven by bias-corrected accelerometer input
- Measurement: `[gps_x, gps_y, baro_z]`

This is not intended to be a production navigation model. It is a shared, deterministic fixture whose job is to answer one question cleanly: does the C EKF match a trusted Python-library reference when both are given the same model, tuning, and sensor sequence?

## Preprocessing and Assumptions

- The first valid GPS fix defines the local origin.
- Latitude/longitude are converted to local meters with an equirectangular approximation around that origin.
- The first barometer pressure sample defines the zero-altitude reference.
- The first 32 IMU samples are averaged to estimate the accelerometer bias/gravity baseline.
- Measurement updates are applied when new GPS data arrives and barometer data is available.
- Because this dataset is mostly stationary bench data, the expected estimate is a near-constant position with velocities staying close to zero.

## Files

- `test/ekfBenchDataExportTest.c`: Runs the C EKF on the bench dataset and exports a CSV trace.
- `scripts/compare_filterpy_bench.py`: Replays the same dataset through FilterPy and compares the results to the C export.

## Output

The C export writes rows with:

- `step`
- `tick_ms`
- `x0..x5` for state
- `P0..P5` for covariance diagonal

The Python script checks:

- state trace agreement
- covariance-diagonal agreement
- exact step/timestamp alignment

Default tolerances are strict (`1e-6`) because both implementations are intentionally matched to the same predict/update ordering and Joseph-form covariance update.

## Notes

- The comparison script uses FilterPy as the reference library surface, but it performs the update math explicitly so it matches this repository's EKF implementation details exactly.
- This fixture is meant for implementation validation first. If we later want a better physical model, we should keep this fixture as a regression harness and add a second, more realistic benchmark alongside it.
