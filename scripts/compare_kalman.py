#!/usr/bin/env python3
import argparse
import csv

import numpy as np
from filterpy.kalman import KalmanFilter


def run_kalman(measurements, dt):
    kf = KalmanFilter(dim_x=2, dim_z=2)
    kf.x = np.array([[0.0], [1.0]])
    kf.P = np.array([[1.0, 0.0], [0.0, 1.0]])
    kf.F = np.array([[1.0, dt], [0.0, 1.0]])
    kf.H = np.eye(2)
    kf.Q = np.array([[0.01, 0.0], [0.0, 0.01]])
    kf.R = np.array([[0.04, 0.0], [0.0, 0.04]])

    rows = []
    for k, z in enumerate(measurements):
        # Predict
        kf.x = kf.F @ kf.x
        kf.P = kf.F @ kf.P @ kf.F.T + kf.Q

        # Update (match C kalmanFilter.c: P = (I - K H) P)
        z = np.array(z).reshape(2, 1)
        S = kf.H @ kf.P @ kf.H.T + kf.R
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
            raise AssertionError("Failed to invert S in python compare.")

        K = kf.P @ kf.H.T @ inv
        y = z - (kf.H @ kf.x)
        kf.x = kf.x + K @ y
        kf.P = (np.eye(kf.P.shape[0]) - K @ kf.H) @ kf.P

        rows.append(
            [
                k,
                kf.x[0, 0],
                kf.x[1, 0],
                kf.P[0, 0],
                kf.P[0, 1],
                kf.P[1, 0],
                kf.P[1, 1],
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
    parser.add_argument("--csv", required=True, help="CSV file from C export")
    parser.add_argument("--tol", type=float, default=None, help="single tolerance")
    parser.add_argument("--tols", default="1e-6,1e-5,1e-4", help="comma-separated tolerances")
    args = parser.parse_args()

    # Measurements must match C test/export
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

    py_rows = run_kalman(measurements, dt=0.1)
    c_rows = read_csv(args.csv)

    if args.tol is not None:
        tolerances = [args.tol]
    else:
        tolerances = [float(t) for t in args.tols.split(",") if t.strip() != ""]

    for tol in tolerances:
        compare_rows(c_rows, py_rows, tol)
        print(f"Kalman comparison passed (tol={tol}).")


if __name__ == "__main__":
    main()
