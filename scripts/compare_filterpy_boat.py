#!/usr/bin/env python3
import argparse
import csv
import math
import sys

import numpy as np
from filterpy.kalman import ExtendedKalmanFilter


WMM_NMAX = 12


def wmm_load(path):
    g = np.zeros((WMM_NMAX + 1, WMM_NMAX + 1))
    h = np.zeros((WMM_NMAX + 1, WMM_NMAX + 1))
    with open(path, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) != 6:
                continue
            n, m = int(parts[0]), int(parts[1])
            if n <= WMM_NMAX and m <= WMM_NMAX:
                # Match C quantization (scaled to 0.01)
                g[n, m] = round(float(parts[2]) * 100.0) / 100.0
                h[n, m] = round(float(parts[3]) * 100.0) / 100.0
    return g, h


def geodetic_to_ecef(lat, lon, h):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    x = (N + h) * cos_lat * cos_lon
    y = (N + h) * cos_lat * sin_lon
    z = (N * (1.0 - e2) + h) * sin_lat
    return x, y, z


def ecef_to_geodetic(x, y, z):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    b = a * (1.0 - f)
    ep2 = (a * a - b * b) / (b * b)
    p = math.sqrt(x * x + y * y)
    if p == 0.0:
        return None
    lon = math.atan2(y, x)
    theta = math.atan2(z * a, p * b)
    sin_t = math.sin(theta)
    cos_t = math.cos(theta)
    lat = math.atan2(z + ep2 * b * sin_t ** 3, p - e2 * a * cos_t ** 3)
    sin_lat = math.sin(lat)
    N = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    h = p / math.cos(lat) - N
    return lat, lon, h


def ned_to_ecef(lat, lon, n, e, d):
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    x = -sin_lat * cos_lon * n - sin_lon * e - cos_lat * cos_lon * d
    y = -sin_lat * sin_lon * n + cos_lon * e - cos_lat * sin_lon * d
    z = cos_lat * n - sin_lat * d
    return x, y, z


def wmm_field_ned(lat, lon, alt, g, h):
    a = 6378137.0
    b = 6356752.314245
    re = 6371200.0
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    rho = (a * a * cos_lat * cos_lat + b * b * sin_lat * sin_lat)
    z = (b * b * sin_lat) / math.sqrt(rho)
    x = (a * a * cos_lat) / math.sqrt(rho)
    r0 = math.sqrt(x * x + z * z)
    st = sin_lat
    ct = cos_lat
    sr = re / (r0 + alt)

    P = np.zeros((WMM_NMAX + 1, WMM_NMAX + 1))
    dP = np.zeros((WMM_NMAX + 1, WMM_NMAX + 1))
    P[0, 0] = 1.0

    for n in range(1, WMM_NMAX + 1):
        for m in range(0, n + 1):
            if n == m:
                P[n, m] = st * P[n - 1, m - 1]
                dP[n, m] = st * dP[n - 1, m - 1] + ct * P[n - 1, m - 1]
            elif n == 1 or m == n - 1:
                P[n, m] = ct * P[n - 1, m]
                dP[n, m] = ct * dP[n - 1, m] - st * P[n - 1, m]
            else:
                k = ((n - 1) * (n - 1) - m * m) / ((2 * n - 1) * (2 * n - 3))
                P[n, m] = ct * P[n - 1, m] - k * P[n - 2, m]
                dP[n, m] = ct * dP[n - 1, m] - st * P[n - 1, m] - k * dP[n - 2, m]

    Br = 0.0
    Bt = 0.0
    Bp = 0.0
    for n in range(1, WMM_NMAX + 1):
        ar = sr ** (n + 2)
        for m in range(0, n + 1):
            cos_m = math.cos(m * lon)
            sin_m = math.sin(m * lon)
            t = g[n, m] * cos_m + h[n, m] * sin_m
            Br += ar * (n + 1) * t * P[n, m]
            Bt -= ar * t * dP[n, m]
            if m != 0:
                Bp += ar * m * (g[n, m] * sin_m - h[n, m] * cos_m) * P[n, m] / ct

    n_comp = -Bt
    e_comp = Bp
    d_comp = -Br
    return n_comp, e_comp, d_comp


def quat_normalize(q):
    norm = np.linalg.norm(q)
    if norm <= 0.0:
        return np.array([1.0, 0.0, 0.0, 0.0])
    return q / norm


def quat_rotate_body_to_ecef(q, v):
    qw, qx, qy, qz = q
    vx, vy, vz = v
    ix = qw * vx + qy * vz - qz * vy
    iy = qw * vy + qz * vx - qx * vz
    iz = qw * vz + qx * vy - qy * vx
    iw = -qx * vx - qy * vy - qz * vz
    ox = ix * qw + iw * -qx + iy * -qz - iz * -qy
    oy = iy * qw + iw * -qy + iz * -qx - ix * -qz
    oz = iz * qw + iw * -qz + ix * -qy - iy * -qx
    return np.array([ox, oy, oz])


def quat_rotate_ecef_to_body(q, v):
    qw, qx, qy, qz = q
    vx, vy, vz = v
    ix = qw * vx - qy * vz + qz * vy
    iy = qw * vy - qz * vx + qx * vz
    iz = qw * vz - qx * vy + qy * vx
    iw = qx * vx + qy * vy + qz * vz
    ox = ix * qw + iw * qx + iy * qz - iz * qy
    oy = iy * qw + iw * qy + iz * qx - ix * qz
    oz = iz * qw + iw * qz + ix * qy - iy * qx
    return np.array([ox, oy, oz])


def boat_fx(x, u, dt):
    K_throttle = 20.0
    K_delta = 0.6
    drag = 0.05
    mass = 50.0 * 0.45359237

    px, py, pz = x[0, 0], x[1, 0], x[2, 0]
    vx, vy, vz = x[3, 0], x[4, 0], x[5, 0]
    q = x[6:10, 0]
    bg = x[10:13, 0]
    ba = x[13:16, 0]
    bm = x[16:19, 0]
    sm = x[19:22, 0]

    omega = np.array([u["gyro_x"], u["gyro_y"], u["gyro_z"]]) - bg
    omega[2] += u["rudder_angle"] * K_delta

    qw, qx, qy, qz = q
    half_dt = 0.5 * dt
    dq_w = -half_dt * (qx * omega[0] + qy * omega[1] + qz * omega[2])
    dq_x =  half_dt * (qw * omega[0] + qy * omega[2] - qz * omega[1])
    dq_y =  half_dt * (qw * omega[1] - qx * omega[2] + qz * omega[0])
    dq_z =  half_dt * (qw * omega[2] + qx * omega[1] - qy * omega[0])
    q = quat_normalize(np.array([qw + dq_w, qx + dq_x, qy + dq_y, qz + dq_z]))

    thrust = u["motor_power"] * K_throttle
    acc_body = np.array([thrust / mass, 0.0, 0.0])
    acc_ecef = quat_rotate_body_to_ecef(q, acc_body)

    vx += (acc_ecef[0] - drag * vx) * dt
    vy += (acc_ecef[1] - drag * vy) * dt
    vz += (acc_ecef[2] - drag * vz) * dt

    px += vx * dt
    py += vy * dt
    pz += vz * dt

    x_out = np.zeros_like(x)
    x_out[0:6, 0] = [px, py, pz, vx, vy, vz]
    x_out[6:10, 0] = q
    x_out[10:13, 0] = bg
    x_out[13:16, 0] = ba
    x_out[16:19, 0] = bm
    x_out[19:22, 0] = sm
    return x_out


def boat_hx(x, g, h):
    g0 = 9.80665
    px, py, pz = x[0, 0], x[1, 0], x[2, 0]
    vx, vy, vz = x[3, 0], x[4, 0], x[5, 0]
    q = x[6:10, 0]
    ba = x[13:16, 0]
    bm = x[16:19, 0]
    sm = x[19:22, 0]

    z = np.zeros((12, 1))
    z[0:6, 0] = [px, py, pz, vx, vy, vz]

    r = math.sqrt(px * px + py * py + pz * pz)
    gx = (-g0 * px / r) if r > 0.0 else 0.0
    gy = (-g0 * py / r) if r > 0.0 else 0.0
    gz = (-g0 * pz / r) if r > 0.0 else -g0
    acc_ecef = np.array([0.0 - gx, 0.0 - gy, 0.0 - gz])
    acc_body = quat_rotate_ecef_to_body(q, acc_ecef)
    z[6:9, 0] = acc_body + ba

    geo = ecef_to_geodetic(px, py, pz)
    if geo:
        lat, lon, alt = geo
        Bn, Be, Bd = wmm_field_ned(lat, lon, alt, g, h)
        Bx, By, Bz = ned_to_ecef(lat, lon, Bn, Be, Bd)
        mag_body = quat_rotate_ecef_to_body(q, np.array([Bx, By, Bz]))
        z[9:12, 0] = sm * mag_body + bm
    return z


def boat_hx_no_wmm(x):
    g0 = 9.80665
    px, py, pz = x[0, 0], x[1, 0], x[2, 0]
    vx, vy, vz = x[3, 0], x[4, 0], x[5, 0]
    q = x[6:10, 0]
    ba = x[13:16, 0]
    bm = x[16:19, 0]
    sm = x[19:22, 0]

    z = np.zeros((12, 1))
    z[0:6, 0] = [px, py, pz, vx, vy, vz]

    r = math.sqrt(px * px + py * py + pz * pz)
    gx = (-g0 * px / r) if r > 0.0 else 0.0
    gy = (-g0 * py / r) if r > 0.0 else 0.0
    gz = (-g0 * pz / r) if r > 0.0 else -g0
    acc_ecef = np.array([0.0 - gx, 0.0 - gy, 0.0 - gz])
    acc_body = quat_rotate_ecef_to_body(q, acc_ecef)
    z[6:9, 0] = acc_body + ba

    geo = ecef_to_geodetic(px, py, pz)
    if geo:
        lat, lon, _ = geo
        Bn, Be, Bd = 0.2, 0.0, 0.45
        Bx, By, Bz = ned_to_ecef(lat, lon, Bn, Be, Bd)
        mag_body = quat_rotate_ecef_to_body(q, np.array([Bx, By, Bz]))
        z[9:12, 0] = sm * mag_body + bm
    return z


def numerical_jacobian(func, x, eps, *args):
    n = x.shape[0]
    y0 = func(x, *args)
    m = y0.shape[0]
    J = np.zeros((m, n))
    for i in range(n):
        xp = x.copy()
        xp[i, 0] += eps
        yp = func(xp, *args)
        J[:, i] = ((yp - y0) / eps).reshape(m)
    return J


def read_csv(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            state = [float(row[f"x{i}"]) for i in range(22)]
            pdiag = [float(row[f"P{i}"]) for i in range(22)]
            rows.append((int(row["step"]), np.array(state), np.array(pdiag)))
    return rows


def compare_rows(c_rows, py_rows, tol_state, tol_p):
    if len(c_rows) != len(py_rows):
        raise AssertionError("Row count mismatch")
    max_state = 0.0
    max_p = 0.0
    for i in range(len(c_rows)):
        c_step, c_state, c_p = c_rows[i]
        p_step, p_state, p_p = py_rows[i]
        if c_step != p_step:
            raise AssertionError(f"Step mismatch at {i}: {c_step} vs {p_step}")
        state_err = np.max(np.abs(c_state - p_state))
        p_err = np.max(np.abs(c_p - p_p))
        max_state = max(max_state, state_err)
        max_p = max(max_p, p_err)
        if state_err > tol_state:
            raise AssertionError(f"State mismatch at step {c_step}")
        if p_err > tol_p:
            raise AssertionError(f"P diag mismatch at step {c_step}")
    return max_state, max_p


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="CSV file from C export")
    parser.add_argument("--wmm", default="WMM.COF", help="WMM.COF path")
    parser.add_argument("--no-wmm", action="store_true", help="disable WMM and use fixed field")
    parser.add_argument("--tol", type=float, default=None, help="single tolerance")
    parser.add_argument("--tols", default="1e-3,1e-2,1e-1", help="comma-separated tolerances")
    parser.add_argument("--tol-state", type=float, default=None, help="state tolerance override")
    parser.add_argument("--tol-p", type=float, default=None, help="P diag tolerance override")
    args = parser.parse_args()

    g = h = None
    if not args.no_wmm:
        g, h = wmm_load(args.wmm)

    ekf = ExtendedKalmanFilter(dim_x=22, dim_z=12)
    x0 = np.zeros((22, 1))
    lat0 = 47.6205 * math.pi / 180.0
    lon0 = -122.3493 * math.pi / 180.0
    x0_e, y0_e, z0_e = geodetic_to_ecef(lat0, lon0, 5.0)
    x0[0, 0] = x0_e
    x0[1, 0] = y0_e
    x0[2, 0] = z0_e
    x0[6, 0] = 1.0
    x0[19, 0] = 1.0
    x0[20, 0] = 1.0
    x0[21, 0] = 1.0
    ekf.x = x0
    ekf.P = np.eye(22)
    ekf.Q = np.eye(22) * 1e-4
    for i in range(10, 22):
        ekf.Q[i, i] = 1e-6
    ekf.R = np.zeros((12, 12))
    diag = [4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.01, 0.01, 0.01]
    for i in range(12):
        ekf.R[i, i] = diag[i]

    c_rows = read_csv(args.csv)

    rows = []
    u = {"gyro_x": 0.0, "gyro_y": 0.0, "gyro_z": 0.01, "motor_power": 0.6, "rudder_angle": 0.05}
    dt = 0.1

    for step in range(len(c_rows)):
        ekf.x = boat_fx(ekf.x, u, dt)
        F = numerical_jacobian(lambda xx: boat_fx(xx, u, dt), ekf.x, 1e-4)
        ekf.P = F @ ekf.P @ F.T + ekf.Q

        if args.no_wmm:
            z = boat_hx_no_wmm(ekf.x)
        else:
            z = boat_hx(ekf.x, g, h)
        noise_gps = 0.5 * math.sin(0.1 * step)
        noise_accel = 0.02 * math.cos(0.2 * step)
        noise_mag = 0.001 * math.sin(0.3 * step)
        gps_drop = (step % 10) < 3
        gps_spike = 50.0 if (step % 15 == 0) else 0.0

        z[0, 0] += noise_gps
        z[1, 0] -= noise_gps
        z[2, 0] += noise_gps
        z[3, 0] += noise_gps
        z[4, 0] -= noise_gps
        z[5, 0] += noise_gps
        z[6, 0] += noise_accel
        z[7, 0] -= noise_accel
        z[8, 0] += noise_accel
        z[9, 0] += noise_mag
        z[10, 0] -= noise_mag
        z[11, 0] += noise_mag

        if gps_drop:
            z[0:6, 0] = ekf.x[0:6, 0]
        elif gps_spike != 0.0:
            z[0, 0] += gps_spike
            z[1, 0] -= gps_spike
            z[2, 0] += gps_spike

        if args.no_wmm:
            H = numerical_jacobian(lambda xx: boat_hx_no_wmm(xx), ekf.x, 1e-4)
            hx = boat_hx_no_wmm(ekf.x)
        else:
            H = numerical_jacobian(lambda xx: boat_hx(xx, g, h), ekf.x, 1e-4)
            hx = boat_hx(ekf.x, g, h)
        y = z - hx
        S = H @ ekf.P @ H.T + ekf.R
        # Jittered inversion to match C behavior
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
        K = ekf.P @ H.T @ inv
        ekf.x = ekf.x + K @ y
        I = np.eye(ekf.P.shape[0])
        ekf.P = (I - K @ H) @ ekf.P @ (I - K @ H).T + K @ ekf.R @ K.T

        rows.append((step, ekf.x[:, 0].copy(), np.diag(ekf.P).copy()))

    if args.tol_state is not None or args.tol_p is not None:
        tol_state = args.tol_state if args.tol_state is not None else 1e-2
        tol_p = args.tol_p if args.tol_p is not None else 2e-2
        max_state, max_p = compare_rows(c_rows, rows, tol_state, tol_p)
        print(f"Filterpy boat comparison passed (tol_state={tol_state}, tol_p={tol_p}). max_state={max_state:.6g} max_p={max_p:.6g}")
        return

    if args.tol is not None:
        tolerances = [args.tol]
    else:
        tolerances = [float(t) for t in args.tols.split(",") if t.strip() != ""]

    any_pass = False
    last_err = None
    for tol in tolerances:
        try:
            max_state, max_p = compare_rows(c_rows, rows, tol, tol)
            print(f"Filterpy boat comparison passed (tol={tol}). max_state={max_state:.6g} max_p={max_p:.6g}")
            any_pass = True
        except AssertionError as err:
            print(f"Filterpy boat comparison failed (tol={tol}): {err}")
            last_err = err
    if not any_pass:
        raise AssertionError(last_err)


if __name__ == "__main__":
    main()
