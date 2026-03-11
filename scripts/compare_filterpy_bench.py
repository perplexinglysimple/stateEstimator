#!/usr/bin/env python3
import argparse
import csv
import math
import sys

import numpy as np


STATE_DIM = 6
MEAS_DIM = 3
CALIBRATION_SAMPLES = 32


def read_imu_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "tick_ms": int(row["tick_ms"]),
                    "ax": float(row["accel_x_mps2"]),
                    "ay": float(row["accel_y_mps2"]),
                    "az": float(row["accel_z_mps2"]),
                }
            )
    return rows


def read_baro_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "tick_ms": int(row["tick_ms"]),
                    "pressure_pa": float(row["pressure_pa"]),
                }
            )
    return rows


def read_gps_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "tick_ms": int(row["tick_ms"]),
                    "valid_fix": int(row["valid_fix"]),
                    "latitude_deg": float(row["latitude_deg"]),
                    "longitude_deg": float(row["longitude_deg"]),
                    "altitude_m": float(row["altitude_m"]),
                }
            )
    return rows


def read_c_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            state = np.array([float(row[f"x{i}"]) for i in range(STATE_DIM)])
            pdiag = np.array([float(row[f"P{i}"]) for i in range(STATE_DIM)])
            rows.append((int(row["step"]), int(row["tick_ms"]), state, pdiag))
    return rows


def geodetic_to_local_meters(lat_deg, lon_deg, lat0_deg, lon0_deg):
    earth_radius_m = 6378137.0
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    lat0 = math.radians(lat0_deg)
    lon0 = math.radians(lon0_deg)
    east_m = earth_radius_m * math.cos(lat0) * (lon - lon0)
    north_m = earth_radius_m * (lat - lat0)
    return east_m, north_m


def pressure_to_relative_altitude(pressure_pa, pressure0_pa):
    if pressure_pa <= 0.0 or pressure0_pa <= 0.0:
        return 0.0
    return 44330.0 * (1.0 - (pressure_pa / pressure0_pa) ** 0.19029495718363465)


def measurement_function(x):
    return x[:MEAS_DIM].copy()


def measurement_jacobian(_x):
    H = np.zeros((MEAS_DIM, STATE_DIM))
    H[0, 0] = 1.0
    H[1, 1] = 1.0
    H[2, 2] = 1.0
    return H


def transition_function(x, ax, ay, az, dt):
    xp = x.copy()
    xp[0, 0] = x[0, 0] + x[3, 0] * dt + 0.5 * ax * dt * dt
    xp[1, 0] = x[1, 0] + x[4, 0] * dt + 0.5 * ay * dt * dt
    xp[2, 0] = x[2, 0] + x[5, 0] * dt + 0.5 * az * dt * dt
    xp[3, 0] = x[3, 0] + ax * dt
    xp[4, 0] = x[4, 0] + ay * dt
    xp[5, 0] = x[5, 0] + az * dt
    return xp


def transition_jacobian(dt):
    F = np.eye(STATE_DIM)
    F[0, 3] = dt
    F[1, 4] = dt
    F[2, 5] = dt
    return F


def compare_rows(c_rows, py_rows, tol_state, tol_p):
    if len(c_rows) != len(py_rows):
        raise AssertionError(f"Row count mismatch: {len(c_rows)} vs {len(py_rows)}")
    max_state = 0.0
    max_p = 0.0
    for c_row, py_row in zip(c_rows, py_rows):
        c_step, c_tick, c_state, c_p = c_row
        py_step, py_tick, py_state, py_p = py_row
        if c_step != py_step or c_tick != py_tick:
            raise AssertionError(
                f"Index mismatch: C(step={c_step}, tick={c_tick}) vs Py(step={py_step}, tick={py_tick})"
            )
        state_err = float(np.max(np.abs(c_state - py_state)))
        p_err = float(np.max(np.abs(c_p - py_p)))
        max_state = max(max_state, state_err)
        max_p = max(max_p, p_err)
        if state_err > tol_state:
            raise AssertionError(f"State mismatch at step {c_step}: {state_err} > {tol_state}")
        if p_err > tol_p:
            raise AssertionError(f"Covariance mismatch at step {c_step}: {p_err} > {tol_p}")
    return max_state, max_p


def run_reference(imu_rows, baro_rows, gps_rows):
    try:
        from filterpy.kalman import ExtendedKalmanFilter
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "filterpy is required for this comparison. Install it with: "
            "python3 -m pip install numpy filterpy"
        ) from exc

    ekf = ExtendedKalmanFilter(dim_x=STATE_DIM, dim_z=MEAS_DIM)
    ekf.x = np.zeros((STATE_DIM, 1))
    ekf.P = np.diag([25.0, 25.0, 25.0, 4.0, 4.0, 4.0])
    ekf.Q = np.diag([0.05, 0.05, 0.05, 0.2, 0.2, 0.2])
    ekf.R = np.diag([9.0, 9.0, 4.0])

    accel_bias = np.zeros(3)
    bias_samples = 0
    previous_tick = None

    pressure0 = None
    latest_baro_z = 0.0
    have_baro = False

    origin_lat = None
    origin_lon = None
    latest_gps_x = 0.0
    latest_gps_y = 0.0
    have_gps = False

    baro_index = 0
    gps_index = 0
    rows = []
    step = 0

    for imu_row in imu_rows:
        measurement_dirty = False
        while baro_index < len(baro_rows) and baro_rows[baro_index]["tick_ms"] <= imu_row["tick_ms"]:
            baro_row = baro_rows[baro_index]
            if pressure0 is None:
                pressure0 = baro_row["pressure_pa"]
            latest_baro_z = pressure_to_relative_altitude(baro_row["pressure_pa"], pressure0)
            have_baro = True
            measurement_dirty = True
            baro_index += 1

        while gps_index < len(gps_rows) and gps_rows[gps_index]["tick_ms"] <= imu_row["tick_ms"]:
            gps_row = gps_rows[gps_index]
            if (
                gps_row["valid_fix"]
                and gps_row["latitude_deg"] != 0.0
                and gps_row["longitude_deg"] != 0.0
            ):
                if origin_lat is None:
                    origin_lat = gps_row["latitude_deg"]
                    origin_lon = gps_row["longitude_deg"]
                latest_gps_x, latest_gps_y = geodetic_to_local_meters(
                    gps_row["latitude_deg"], gps_row["longitude_deg"], origin_lat, origin_lon
                )
                have_gps = True
                measurement_dirty = True
            gps_index += 1

        if bias_samples < CALIBRATION_SAMPLES:
            accel_bias += np.array([imu_row["ax"], imu_row["ay"], imu_row["az"]])
            bias_samples += 1

        if previous_tick is None:
            previous_tick = imu_row["tick_ms"]
            continue

        dt = (imu_row["tick_ms"] - previous_tick) / 1000.0
        previous_tick = imu_row["tick_ms"]
        if dt <= 0.0:
            continue

        accel = np.array([[imu_row["ax"]], [imu_row["ay"]], [imu_row["az"]]])
        accel_bias_mean = accel_bias / max(bias_samples, 1)
        dynamic_accel = accel.reshape(3) - accel_bias_mean

        ekf.x = transition_function(ekf.x, dynamic_accel[0], dynamic_accel[1], dynamic_accel[2], dt)
        ekf.F = transition_jacobian(dt)
        ekf.P = ekf.F @ ekf.P @ ekf.F.T + ekf.Q

        if measurement_dirty and have_baro and have_gps:
            H = measurement_jacobian(ekf.x)
            z = np.array([[latest_gps_x], [latest_gps_y], [latest_baro_z]])
            hx = measurement_function(ekf.x)
            y = z - hx
            S = H @ ekf.P @ H.T + ekf.R
            inv = None
            jitter = 1e-6
            for _ in range(3):
                try:
                    inv = np.linalg.inv(S)
                    break
                except np.linalg.LinAlgError:
                    S = S + jitter * np.eye(S.shape[0])
                    jitter *= 100.0
            if inv is None:
                raise AssertionError("Failed to invert innovation matrix S in python reference.")
            K = ekf.P @ H.T @ inv
            ekf.x = ekf.x + K @ y
            I = np.eye(STATE_DIM)
            ekf.P = (I - K @ H) @ ekf.P @ (I - K @ H).T + K @ ekf.R @ K.T

        rows.append((step, imu_row["tick_ms"], ekf.x[:, 0].copy(), np.diag(ekf.P).copy()))
        step += 1

    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="CSV exported by ekfBenchDataExportTest")
    parser.add_argument("--imu", default="test_data/ICM20948.CSV", help="IMU CSV path")
    parser.add_argument("--gps", default="test_data/SAM_M8Q.CSV", help="GPS CSV path")
    parser.add_argument("--baro", default="test_data/BMP384.CSV", help="Barometer CSV path")
    parser.add_argument("--tol", type=float, default=None, help="Single tolerance for state and covariance")
    parser.add_argument("--tol-state", type=float, default=1e-6, help="State comparison tolerance")
    parser.add_argument("--tol-p", type=float, default=1e-6, help="Covariance comparison tolerance")
    args = parser.parse_args()

    imu_rows = read_imu_rows(args.imu)
    baro_rows = read_baro_rows(args.baro)
    gps_rows = read_gps_rows(args.gps)
    c_rows = read_c_rows(args.csv)
    py_rows = run_reference(imu_rows, baro_rows, gps_rows)

    tol_state = args.tol if args.tol is not None else args.tol_state
    tol_p = args.tol if args.tol is not None else args.tol_p
    max_state, max_p = compare_rows(c_rows, py_rows, tol_state, tol_p)
    print(
        f"FilterPy bench comparison passed (tol_state={tol_state}, tol_p={tol_p}). "
        f"max_state={max_state:.6g} max_p={max_p:.6g}"
    )


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        raise
