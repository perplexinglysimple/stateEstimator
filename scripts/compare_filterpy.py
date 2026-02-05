#!/usr/bin/env python3
import argparse
import csv
import math
import sys

import numpy as np
from filterpy.kalman import ExtendedKalmanFilter



def fx_linear(x, dt):
    F = np.array([[1.0, dt], [0.0, 1.0]])
    return F @ x


def F_linear(x, dt):
    return np.array([[1.0, dt], [0.0, 1.0]])


def hx_linear(x):
    return x.copy()


def H_linear(x):
    return np.eye(2)


def fx_nonlinear(x, dt):
    x0 = x[0, 0]
    x1 = x[1, 0]
    y0 = x0 + x1 * dt + 0.5 * dt * dt * math.sin(x0)
    y1 = x1 + dt * math.cos(x1)
    return np.array([[y0], [y1]])


def F_nonlinear(x, dt):
    x0 = x[0, 0]
    x1 = x[1, 0]
    return np.array(
        [
            [1.0 + 0.5 * dt * dt * math.cos(x0), dt],
            [0.0, 1.0 - dt * math.sin(x1)],
        ]
    )


def hx_nonlinear(x):
    x0 = x[0, 0]
    x1 = x[1, 0]
    return np.array([[x0 * x0], [x1 * x1]])


def H_nonlinear(x):
    x0 = x[0, 0]
    x1 = x[1, 0]
    return np.array([[2.0 * x0, 0.0], [0.0, 2.0 * x1]])


def run_filterpy(measurements, dt, model):
    ekf = ExtendedKalmanFilter(dim_x=2, dim_z=2)
    ekf.x = np.array([[0.0], [1.0]])
    ekf.P = np.array([[1.0, 0.0], [0.0, 1.0]])
    ekf.Q = np.array([[0.01, 0.0], [0.0, 0.01]])
    ekf.R = np.array([[0.04, 0.0], [0.0, 0.04]])

    rows = []
    for k, z in enumerate(measurements):
        if model == "linear":
            ekf.F = F_linear(ekf.x, dt)
            ekf.predict()
            ekf.update(np.array(z).reshape(2, 1), H_linear, hx_linear)
        else:
            # Match C EKF behavior: propagate with f(x), then use F at predicted state
            ekf.x = fx_nonlinear(ekf.x, dt)
            ekf.F = F_nonlinear(ekf.x, dt)
            ekf.P = ekf.F @ ekf.P @ ekf.F.T + ekf.Q
            ekf.update(np.array(z).reshape(2, 1), H_nonlinear, hx_nonlinear)
        rows.append(
            [
                k,
                ekf.x[0, 0],
                ekf.x[1, 0],
                ekf.P[0, 0],
                ekf.P[0, 1],
                ekf.P[1, 0],
                ekf.P[1, 1],
            ]
        )
    return rows


def read_csv(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                [
                    int(row["step"]),
                    float(row["x0"]),
                    float(row["x1"]),
                    float(row["P00"]),
                    float(row["P01"]),
                    float(row["P10"]),
                    float(row["P11"]),
                ]
            )
    return rows


def compare_rows(c_rows, py_rows, tol):
    if len(c_rows) != len(py_rows):
        raise AssertionError("Row count mismatch")
    for i in range(len(c_rows)):
        c = c_rows[i]
        p = py_rows[i]
        if c[0] != p[0]:
            raise AssertionError(f"Step mismatch at {i}: {c[0]} vs {p[0]}")
        for j in range(1, 7):
            if abs(c[j] - p[j]) > tol:
                raise AssertionError(
                    f"Mismatch at step {c[0]} index {j}: c={c[j]} py={p[j]}"
                )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="CSV file from C test")
    parser.add_argument("--tol", type=float, default=None, help="single tolerance")
    parser.add_argument("--tols", default="1e-6,1e-5,1e-4", help="comma-separated tolerances")
    parser.add_argument("--model", choices=["linear", "nonlinear"], default="linear")
    args = parser.parse_args()

    # Measurements must match C test
    if args.model == "linear":
        measurements = [
            [0.0, 1.0],
            [0.1, 1.0],
            [0.2, 1.0],
            [0.3, 1.0],
            [0.4, 1.0],
            [0.5, 1.0],
            [0.6, 1.0],
            [0.7, 1.0],
            [0.8, 1.0],
            [0.9, 1.0],
        ]
    else:
        measurements = [
            [0.00, 1.00],
            [0.05, 0.98],
            [0.10, 0.97],
            [0.16, 0.95],
            [0.21, 0.94],
            [0.27, 0.92],
            [0.33, 0.90],
            [0.39, 0.89],
            [0.46, 0.87],
            [0.52, 0.86],
        ]

    py_rows = run_filterpy(measurements, dt=0.1, model=args.model)
    c_rows = read_csv(args.csv)

    if args.tol is not None:
        tolerances = [args.tol]
    else:
        tolerances = [float(t) for t in args.tols.split(",") if t.strip() != ""]

    for tol in tolerances:
        compare_rows(c_rows, py_rows, tol)
        print(f"Filterpy comparison passed (model={args.model}, tol={tol}).")


if __name__ == "__main__":
    main()
