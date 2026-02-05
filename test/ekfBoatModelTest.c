#include "EKF.h"
#include <math.h>
#include <stdio.h>

typedef struct BoatInputs_
{
    double gyro_x;
    double gyro_y;
    double gyro_z;
    double motor_power;  // 0..1
    double rudder_angle; // radians
} BoatInputs;

static void BoatTransition(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void BoatMeasurement(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void BoatUpdateA(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

static void geodeticToECEF(double lat_rad, double lon_rad, double h_m, double* x, double* y, double* z);
static int ecefToGeodetic(double x, double y, double z, double* lat_rad, double* lon_rad, double* h_m);
static void nedToEcef(double lat_rad, double lon_rad, double n, double e, double d, double* x, double* y, double* z);
static int wmmLoadModel(const char* path);
static int wmmFieldNED(double lat_rad, double lon_rad, double alt_m, double* n, double* e, double* d);
static void quatNormalize(double* qw, double* qx, double* qy, double* qz);
static void quatRotateBodyToECEF(double qw, double qx, double qy, double qz, double vx, double vy, double vz,
                                 double* ox, double* oy, double* oz);
static void quatRotateECEFToBody(double qw, double qx, double qy, double qz, double vx, double vy, double vz,
                                 double* ox, double* oy, double* oz);

#define WMM_NMAX 12
static int g_coeff[WMM_NMAX + 1][WMM_NMAX + 1];
static int h_coeff[WMM_NMAX + 1][WMM_NMAX + 1];
static int wmm_loaded = 0;

#define BOAT_STATE_DIM 22
#define BOAT_MEAS_DIM 12

int main()
{
    const int n = BOAT_STATE_DIM;
    const int m = BOAT_MEAS_DIM;
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    INIT_MATRIX(ekf.x, BOAT_STATE_DIM, 1);
    INIT_MATRIX(ekf.P, BOAT_STATE_DIM, BOAT_STATE_DIM);
    INIT_MATRIX(ekf.Q, BOAT_STATE_DIM, BOAT_STATE_DIM);
    INIT_MATRIX(ekf.R, BOAT_MEAS_DIM, BOAT_MEAS_DIM);
    INIT_MATRIX(ekf.A, BOAT_STATE_DIM, BOAT_STATE_DIM);

    STATIC_MATRIX_DIRECTIVE(options.x0, BOAT_STATE_DIM, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, BOAT_STATE_DIM, BOAT_STATE_DIM, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, BOAT_STATE_DIM, BOAT_STATE_DIM, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, BOAT_MEAS_DIM, BOAT_MEAS_DIM, R);
    STATIC_MATRIX_DIRECTIVE(options.A, BOAT_STATE_DIM, BOAT_STATE_DIM, A);

    // Initial state from geodetic
    const double kPi = 3.14159265358979323846;
    double lat0 = 47.6205 * kPi / 180.0;
    double lon0 = -122.3493 * kPi / 180.0;
    double h0 = 5.0;
    double x0_e, y0_e, z0_e;
    geodeticToECEF(lat0, lon0, h0, &x0_e, &y0_e, &z0_e);
    SET_MATRIX(*(options.x0), 0, 0, x0_e);
    SET_MATRIX(*(options.x0), 1, 0, y0_e);
    SET_MATRIX(*(options.x0), 2, 0, z0_e);
    SET_MATRIX(*(options.x0), 3, 0, 0.0);
    SET_MATRIX(*(options.x0), 4, 0, 0.0);
    SET_MATRIX(*(options.x0), 5, 0, 0.0);
    // Quaternion identity
    SET_MATRIX(*(options.x0), 6, 0, 1.0);
    SET_MATRIX(*(options.x0), 7, 0, 0.0);
    SET_MATRIX(*(options.x0), 8, 0, 0.0);
    SET_MATRIX(*(options.x0), 9, 0, 0.0);
    // gyro bias (10..12), accel bias (13..15), mag bias (16..18), mag scale (19..21)
    for (int i = 10; i < BOAT_STATE_DIM; ++i)
    {
        SET_MATRIX(*(options.x0), i, 0, 0.0);
    }
    SET_MATRIX(*(options.x0), 19, 0, 1.0);
    SET_MATRIX(*(options.x0), 20, 0, 1.0);
    SET_MATRIX(*(options.x0), 21, 0, 1.0);

    for (int i = 0; i < n; ++i)
    {
        SET_MATRIX(*(options.P0), i, i, 1.0);
        SET_MATRIX(*(options.Q), i, i, 1e-4);
    }
    // Larger uncertainty for bias/random-walk states
    for (int i = 10; i < BOAT_STATE_DIM; ++i)
    {
        SET_MATRIX(*(options.Q), i, i, 1e-6);
    }
    // Measurement noise: GPS pos/vel, accel, mag
    SET_MATRIX(*(options.R), 0, 0, 4.0);
    SET_MATRIX(*(options.R), 1, 1, 4.0);
    SET_MATRIX(*(options.R), 2, 2, 4.0);
    SET_MATRIX(*(options.R), 3, 3, 1.0);
    SET_MATRIX(*(options.R), 4, 4, 1.0);
    SET_MATRIX(*(options.R), 5, 5, 1.0);
    SET_MATRIX(*(options.R), 6, 6, 0.5);
    SET_MATRIX(*(options.R), 7, 7, 0.5);
    SET_MATRIX(*(options.R), 8, 8, 0.5);
    SET_MATRIX(*(options.R), 9, 9, 0.01);
    SET_MATRIX(*(options.R), 10, 10, 0.01);
    SET_MATRIX(*(options.R), 11, 11, 0.01);
    for (int i = 0; i < n; ++i)
    {
        SET_MATRIX(*(options.A), i, i, 1.0);
    }

    options.n = n;
    options.f = BoatTransition;
    options.h = BoatMeasurement;
    options.updateAMatrix = BoatUpdateA;
    options.numberOfStates = n;
    options.numberOfMeasurements = m;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    if (wmmLoadModel("WMM.COF") != 0 && wmmLoadModel("..\\WMM.COF") != 0)
    {
        LOG_ERROR("Failed to load WMM.COF.");
        return -1;
    }

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in boat model test.");
        return -1;
    }

    BoatInputs inputs = {0};
    inputs.gyro_x = 0.0;
    inputs.gyro_y = 0.0;
    inputs.gyro_z = 0.01;
    inputs.motor_power = 0.6;
    inputs.rudder_angle = 0.05;

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, BOAT_MEAS_DIM, 1, z);

    const double dt = 0.1;
    for (int step = 0; step < 25; ++step)
    {
        if (EKFPredict(&ekf, dt, &inputs) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFPredict() failed in boat model test.");
            return -1;
        }

        // Build measurement using the model (perfect measurement)
        BoatMeasurement(ekf._x_predicted, measurement.z, &ekf, NULL);

        // Inject deterministic noise into measurements
        double noise_gps = 0.5 * sin(0.1 * step);
        double noise_accel = 0.02 * cos(0.2 * step);
        double noise_mag = 0.001 * sin(0.3 * step);

        // Periodic GPS dropout: 3 steps every 10
        bool gps_drop = (step % 10) < 3;
        // Periodic outlier spike
        double gps_spike = (step % 15 == 0) ? 50.0 : 0.0;

        SET_MATRIX(*(measurement.z), 0, 0, ACCESS_MATRIX(*(measurement.z), 0, 0) + noise_gps);
        SET_MATRIX(*(measurement.z), 1, 0, ACCESS_MATRIX(*(measurement.z), 1, 0) - noise_gps);
        SET_MATRIX(*(measurement.z), 2, 0, ACCESS_MATRIX(*(measurement.z), 2, 0) + noise_gps);

        SET_MATRIX(*(measurement.z), 3, 0, ACCESS_MATRIX(*(measurement.z), 3, 0) + noise_gps);
        SET_MATRIX(*(measurement.z), 4, 0, ACCESS_MATRIX(*(measurement.z), 4, 0) - noise_gps);
        SET_MATRIX(*(measurement.z), 5, 0, ACCESS_MATRIX(*(measurement.z), 5, 0) + noise_gps);

        SET_MATRIX(*(measurement.z), 6, 0, ACCESS_MATRIX(*(measurement.z), 6, 0) + noise_accel);
        SET_MATRIX(*(measurement.z), 7, 0, ACCESS_MATRIX(*(measurement.z), 7, 0) - noise_accel);
        SET_MATRIX(*(measurement.z), 8, 0, ACCESS_MATRIX(*(measurement.z), 8, 0) + noise_accel);

        SET_MATRIX(*(measurement.z), 9, 0, ACCESS_MATRIX(*(measurement.z), 9, 0) + noise_mag);
        SET_MATRIX(*(measurement.z), 10, 0, ACCESS_MATRIX(*(measurement.z), 10, 0) - noise_mag);
        SET_MATRIX(*(measurement.z), 11, 0, ACCESS_MATRIX(*(measurement.z), 11, 0) + noise_mag);

        if (gps_drop)
        {
            // Replace GPS pos/vel with predicted values
            SET_MATRIX(*(measurement.z), 0, 0, ACCESS_MATRIX(*(ekf._x_predicted), 0, 0));
            SET_MATRIX(*(measurement.z), 1, 0, ACCESS_MATRIX(*(ekf._x_predicted), 1, 0));
            SET_MATRIX(*(measurement.z), 2, 0, ACCESS_MATRIX(*(ekf._x_predicted), 2, 0));
            SET_MATRIX(*(measurement.z), 3, 0, ACCESS_MATRIX(*(ekf._x_predicted), 3, 0));
            SET_MATRIX(*(measurement.z), 4, 0, ACCESS_MATRIX(*(ekf._x_predicted), 4, 0));
            SET_MATRIX(*(measurement.z), 5, 0, ACCESS_MATRIX(*(ekf._x_predicted), 5, 0));
        }
        else if (gps_spike != 0.0)
        {
            SET_MATRIX(*(measurement.z), 0, 0, ACCESS_MATRIX(*(measurement.z), 0, 0) + gps_spike);
            SET_MATRIX(*(measurement.z), 1, 0, ACCESS_MATRIX(*(measurement.z), 1, 0) - gps_spike);
            SET_MATRIX(*(measurement.z), 2, 0, ACCESS_MATRIX(*(measurement.z), 2, 0) + gps_spike);
        }

        if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFUpdate() failed in boat model test.");
            return -1;
        }

        // Sanity check quaternion normalization
        double qw = ACCESS_MATRIX(*(ekf.x), 6, 0);
        double qx = ACCESS_MATRIX(*(ekf.x), 7, 0);
        double qy = ACCESS_MATRIX(*(ekf.x), 8, 0);
        double qz = ACCESS_MATRIX(*(ekf.x), 9, 0);
        double qnorm = sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
        if (fabs(qnorm - 1.0) > 1e-2)
        {
            LOG_ERROR("Quaternion drifted from unit norm.");
            return -1;
        }

        // GPS conversion round-trip check (ECEF -> geodetic)
        double lat_chk, lon_chk, h_chk;
        if (ecefToGeodetic(ACCESS_MATRIX(*(ekf.x), 0, 0), ACCESS_MATRIX(*(ekf.x), 1, 0), ACCESS_MATRIX(*(ekf.x), 2, 0),
                           &lat_chk, &lon_chk, &h_chk) != 0)
        {
            LOG_ERROR("ECEF to geodetic conversion failed.");
            return -1;
        }
    }

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed in boat model test.");
        return -1;
    }

    return 0;
}

static void BoatTransition(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    const BoatInputs* inputs = (const BoatInputs*) userData;
    const double dt = 0.1;
    const double mass = 50.0 * 0.45359237;
    const double K_throttle = 20.0;
    const double K_delta = 0.6; // rad/s per rad
    const double drag = 0.05;

    double px = ACCESS_MATRIX(*x, 0, 0);
    double py = ACCESS_MATRIX(*x, 1, 0);
    double pz = ACCESS_MATRIX(*x, 2, 0);
    double vx = ACCESS_MATRIX(*x, 3, 0);
    double vy = ACCESS_MATRIX(*x, 4, 0);
    double vz = ACCESS_MATRIX(*x, 5, 0);

    double qw = ACCESS_MATRIX(*x, 6, 0);
    double qx = ACCESS_MATRIX(*x, 7, 0);
    double qy = ACCESS_MATRIX(*x, 8, 0);
    double qz = ACCESS_MATRIX(*x, 9, 0);

    double bgx = ACCESS_MATRIX(*x, 10, 0);
    double bgy = ACCESS_MATRIX(*x, 11, 0);
    double bgz = ACCESS_MATRIX(*x, 12, 0);
    double bax = ACCESS_MATRIX(*x, 13, 0);
    double bay = ACCESS_MATRIX(*x, 14, 0);
    double baz = ACCESS_MATRIX(*x, 15, 0);
    double bmx = ACCESS_MATRIX(*x, 16, 0);
    double bmy = ACCESS_MATRIX(*x, 17, 0);
    double bmz = ACCESS_MATRIX(*x, 18, 0);
    double smx = ACCESS_MATRIX(*x, 19, 0);
    double smy = ACCESS_MATRIX(*x, 20, 0);
    double smz = ACCESS_MATRIX(*x, 21, 0);

    double omega_x = (inputs ? inputs->gyro_x : 0.0) - bgx;
    double omega_y = (inputs ? inputs->gyro_y : 0.0) - bgy;
    double omega_z = (inputs ? inputs->gyro_z : 0.0) - bgz;
    omega_z += (inputs ? inputs->rudder_angle : 0.0) * K_delta;

    // Quaternion integration (body rates)
    double half_dt = 0.5 * dt;
    double dq_w = -half_dt * (qx * omega_x + qy * omega_y + qz * omega_z);
    double dq_x = half_dt * (qw * omega_x + qy * omega_z - qz * omega_y);
    double dq_y = half_dt * (qw * omega_y - qx * omega_z + qz * omega_x);
    double dq_z = half_dt * (qw * omega_z + qx * omega_y - qy * omega_x);

    qw += dq_w;
    qx += dq_x;
    qy += dq_y;
    qz += dq_z;
    quatNormalize(&qw, &qx, &qy, &qz);

    // Thrust in body frame (forward x)
    double thrust = (inputs ? inputs->motor_power : 0.0) * K_throttle;
    double ax_body = thrust / mass;
    double ay_body = 0.0;
    double az_body = 0.0;

    // Rotate acceleration into ECEF
    double ax_e, ay_e, az_e;
    quatRotateBodyToECEF(qw, qx, qy, qz, ax_body, ay_body, az_body, &ax_e, &ay_e, &az_e);

    vx += (ax_e - drag * vx) * dt;
    vy += (ay_e - drag * vy) * dt;
    vz += (az_e - drag * vz) * dt;

    px += vx * dt;
    py += vy * dt;
    pz += vz * dt;

    SET_MATRIX(*x_pred, 0, 0, px);
    SET_MATRIX(*x_pred, 1, 0, py);
    SET_MATRIX(*x_pred, 2, 0, pz);
    SET_MATRIX(*x_pred, 3, 0, vx);
    SET_MATRIX(*x_pred, 4, 0, vy);
    SET_MATRIX(*x_pred, 5, 0, vz);
    SET_MATRIX(*x_pred, 6, 0, qw);
    SET_MATRIX(*x_pred, 7, 0, qx);
    SET_MATRIX(*x_pred, 8, 0, qy);
    SET_MATRIX(*x_pred, 9, 0, qz);
    SET_MATRIX(*x_pred, 10, 0, bgx);
    SET_MATRIX(*x_pred, 11, 0, bgy);
    SET_MATRIX(*x_pred, 12, 0, bgz);
    SET_MATRIX(*x_pred, 13, 0, bax);
    SET_MATRIX(*x_pred, 14, 0, bay);
    SET_MATRIX(*x_pred, 15, 0, baz);
    SET_MATRIX(*x_pred, 16, 0, bmx);
    SET_MATRIX(*x_pred, 17, 0, bmy);
    SET_MATRIX(*x_pred, 18, 0, bmz);
    SET_MATRIX(*x_pred, 19, 0, smx);
    SET_MATRIX(*x_pred, 20, 0, smy);
    SET_MATRIX(*x_pred, 21, 0, smz);
}

static void BoatMeasurement(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    const double g = 9.80665;
    double Bn = 0.0;
    double Be = 0.0;
    double Bd = 0.0;

    double px = ACCESS_MATRIX(*x, 0, 0);
    double py = ACCESS_MATRIX(*x, 1, 0);
    double pz = ACCESS_MATRIX(*x, 2, 0);
    double qw = ACCESS_MATRIX(*x, 6, 0);
    double qx = ACCESS_MATRIX(*x, 7, 0);
    double qy = ACCESS_MATRIX(*x, 8, 0);
    double qz = ACCESS_MATRIX(*x, 9, 0);
    double bax = ACCESS_MATRIX(*x, 13, 0);
    double bay = ACCESS_MATRIX(*x, 14, 0);
    double baz = ACCESS_MATRIX(*x, 15, 0);
    double bmx = ACCESS_MATRIX(*x, 16, 0);
    double bmy = ACCESS_MATRIX(*x, 17, 0);
    double bmz = ACCESS_MATRIX(*x, 18, 0);
    double smx = ACCESS_MATRIX(*x, 19, 0);
    double smy = ACCESS_MATRIX(*x, 20, 0);
    double smz = ACCESS_MATRIX(*x, 21, 0);

    // GPS position (ECEF)
    SET_MATRIX(*z, 0, 0, px);
    SET_MATRIX(*z, 1, 0, py);
    SET_MATRIX(*z, 2, 0, pz);
    SET_MATRIX(*z, 3, 0, ACCESS_MATRIX(*x, 3, 0));
    SET_MATRIX(*z, 4, 0, ACCESS_MATRIX(*x, 4, 0));
    SET_MATRIX(*z, 5, 0, ACCESS_MATRIX(*x, 5, 0));

    // Gravity direction in ECEF (approximate radial)
    double r = sqrt(px * px + py * py + pz * pz);
    double gx = (r > 0.0) ? (-g * px / r) : 0.0;
    double gy = (r > 0.0) ? (-g * py / r) : 0.0;
    double gz = (r > 0.0) ? (-g * pz / r) : -g;

    // Expected accel measurement (specific force)
    double ax_e = 0.0;
    double ay_e = 0.0;
    double az_e = 0.0;
    double ax_b, ay_b, az_b;
    quatRotateECEFToBody(qw, qx, qy, qz, ax_e - gx, ay_e - gy, az_e - gz, &ax_b, &ay_b, &az_b);
    SET_MATRIX(*z, 6, 0, ax_b + bax);
    SET_MATRIX(*z, 7, 0, ay_b + bay);
    SET_MATRIX(*z, 8, 0, az_b + baz);

    // Expected mag measurement in body (WMM NED -> ECEF -> body)
    double lat_rad = 0.0;
    double lon_rad = 0.0;
    double alt_m = 0.0;
    if (ecefToGeodetic(px, py, pz, &lat_rad, &lon_rad, &alt_m) == 0)
    {
        if (wmmFieldNED(lat_rad, lon_rad, alt_m, &Bn, &Be, &Bd) == 0)
        {
            double Bx_e, By_e, Bz_e;
            nedToEcef(lat_rad, lon_rad, Bn, Be, Bd, &Bx_e, &By_e, &Bz_e);
            double mx_b, my_b, mz_b;
            quatRotateECEFToBody(qw, qx, qy, qz, Bx_e, By_e, Bz_e, &mx_b, &my_b, &mz_b);
            SET_MATRIX(*z, 9, 0, smx * mx_b + bmx);
            SET_MATRIX(*z, 10, 0, smy * my_b + bmy);
            SET_MATRIX(*z, 11, 0, smz * mz_b + bmz);
        }
    }
}

static void BoatUpdateA(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    (void) x;
    (void) ekf;
    (void) time;
    (void) userData;
    // This model is nonlinear; A is not used directly.
    for (int i = 0; i < A->row; ++i)
    {
        for (int j = 0; j < A->col; ++j)
        {
            SET_MATRIX(*A, i, j, (i == j) ? 1.0 : 0.0);
        }
    }
}

static void geodeticToECEF(double lat_rad, double lon_rad, double h_m, double* x, double* y, double* z)
{
    const double a = 6378137.0;
    const double f = 1.0 / 298.257223563;
    const double e2 = f * (2.0 - f);
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    *x = (N + h_m) * cos_lat * cos_lon;
    *y = (N + h_m) * cos_lat * sin_lon;
    *z = (N * (1.0 - e2) + h_m) * sin_lat;
}

static int ecefToGeodetic(double x, double y, double z, double* lat_rad, double* lon_rad, double* h_m)
{
    const double a = 6378137.0;
    const double f = 1.0 / 298.257223563;
    const double e2 = f * (2.0 - f);
    const double b = a * (1.0 - f);
    const double ep2 = (a * a - b * b) / (b * b);

    double p = sqrt(x * x + y * y);
    if (p == 0.0)
    {
        return -1;
    }
    *lon_rad = atan2(y, x);
    double theta = atan2(z * a, p * b);
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    *lat_rad = atan2(z + ep2 * b * sin_t * sin_t * sin_t, p - e2 * a * cos_t * cos_t * cos_t);
    double sin_lat = sin(*lat_rad);
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    *h_m = p / cos(*lat_rad) - N;
    return 0;
}

static void nedToEcef(double lat_rad, double lon_rad, double n, double e, double d, double* x, double* y, double* z)
{
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);

    *x = -sin_lat * cos_lon * n - sin_lon * e - cos_lat * cos_lon * d;
    *y = -sin_lat * sin_lon * n + cos_lon * e - cos_lat * sin_lon * d;
    *z = cos_lat * n - sin_lat * d;
}

static void quatNormalize(double* qw, double* qx, double* qy, double* qz)
{
    double norm = sqrt((*qw) * (*qw) + (*qx) * (*qx) + (*qy) * (*qy) + (*qz) * (*qz));
    if (norm <= 0.0)
    {
        *qw = 1.0;
        *qx = 0.0;
        *qy = 0.0;
        *qz = 0.0;
        return;
    }
    *qw /= norm;
    *qx /= norm;
    *qy /= norm;
    *qz /= norm;
}

static void quatRotateBodyToECEF(double qw, double qx, double qy, double qz, double vx, double vy, double vz,
                                 double* ox, double* oy, double* oz)
{
    // q * v * q_conj
    double ix = qw * vx + qy * vz - qz * vy;
    double iy = qw * vy + qz * vx - qx * vz;
    double iz = qw * vz + qx * vy - qy * vx;
    double iw = -qx * vx - qy * vy - qz * vz;

    *ox = ix * qw + iw * -qx + iy * -qz - iz * -qy;
    *oy = iy * qw + iw * -qy + iz * -qx - ix * -qz;
    *oz = iz * qw + iw * -qz + ix * -qy - iy * -qx;
}

static void quatRotateECEFToBody(double qw, double qx, double qy, double qz, double vx, double vy, double vz,
                                 double* ox, double* oy, double* oz)
{
    // q_conj * v * q
    double ix = qw * vx - qy * vz + qz * vy;
    double iy = qw * vy - qz * vx + qx * vz;
    double iz = qw * vz - qx * vy + qy * vx;
    double iw = qx * vx + qy * vy + qz * vz;

    *ox = ix * qw + iw * qx + iy * qz - iz * qy;
    *oy = iy * qw + iw * qy + iz * qx - ix * qz;
    *oz = iz * qw + iw * qz + ix * qy - iy * qx;
}

static int wmmLoadModel(const char* path)
{
    FILE* f = NULL;
#ifdef _WIN32
    if (fopen_s(&f, path, "r") != 0)
    {
        f = NULL;
    }
#else
    f = fopen(path, "r");
#endif
    if (!f)
    {
        return -1;
    }
    for (int n = 0; n <= WMM_NMAX; ++n)
    {
        for (int m = 0; m <= WMM_NMAX; ++m)
        {
            g_coeff[n][m] = 0;
            h_coeff[n][m] = 0;
        }
    }
    char line[256];
    if (!fgets(line, sizeof(line), f))
    {
        fclose(f);
        return -1;
    }
    int n, m;
    double g, h, dg, dh;
    while (fgets(line, sizeof(line), f))
    {
#ifdef _WIN32
        if (sscanf_s(line, "%d %d %lf %lf %lf %lf", &n, &m, &g, &h, &dg, &dh) == 6)
#else
        if (sscanf(line, "%d %d %lf %lf %lf %lf", &n, &m, &g, &h, &dg, &dh) == 6)
#endif
        {
            if (n <= WMM_NMAX && m <= WMM_NMAX)
            {
                g_coeff[n][m] = (int) lrint(g * 100.0);
                h_coeff[n][m] = (int) lrint(h * 100.0);
            }
        }
    }
    fclose(f);
    wmm_loaded = 1;
    return 0;
}

static int wmmFieldNED(double lat_rad, double lon_rad, double alt_m, double* n, double* e, double* d)
{
    if (!wmm_loaded)
    {
        return -1;
    }
    const double a = 6378137.0;
    const double b = 6356752.314245;
    const double re = 6371200.0;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);
    double sin_lon = sin(lon_rad);
    double cos_lon = cos(lon_rad);

    double r = sqrt((a * a * cos_lat * cos_lat + b * b * sin_lat * sin_lat));
    (void) sin_lon;
    (void) cos_lon;
    (void) r;
    double rho = (a * a * cos_lat * cos_lat + b * b * sin_lat * sin_lat);
    double z = (b * b * sin_lat) / sqrt(rho);
    double x = (a * a * cos_lat) / sqrt(rho);
    double r0 = sqrt(x * x + z * z);
    double st = sin_lat;
    double ct = cos_lat;
    double sr = (re / (r0 + alt_m));

    double P[WMM_NMAX + 1][WMM_NMAX + 1] = {0};
    double dP[WMM_NMAX + 1][WMM_NMAX + 1] = {0};
    P[0][0] = 1.0;

    for (int n_idx = 1; n_idx <= WMM_NMAX; ++n_idx)
    {
        for (int m = 0; m <= n_idx; ++m)
        {
            if (n_idx == m)
            {
                P[n_idx][m] = st * P[n_idx - 1][m - 1];
                dP[n_idx][m] = st * dP[n_idx - 1][m - 1] + ct * P[n_idx - 1][m - 1];
            }
            else if (n_idx == 1 || m == n_idx - 1)
            {
                P[n_idx][m] = ct * P[n_idx - 1][m];
                dP[n_idx][m] = ct * dP[n_idx - 1][m] - st * P[n_idx - 1][m];
            }
            else
            {
                double k = ((n_idx - 1) * (n_idx - 1) - m * m) / (double) ((2 * n_idx - 1) * (2 * n_idx - 3));
                P[n_idx][m] = ct * P[n_idx - 1][m] - k * P[n_idx - 2][m];
                dP[n_idx][m] = ct * dP[n_idx - 1][m] - st * P[n_idx - 1][m] - k * dP[n_idx - 2][m];
            }
        }
    }

    double Br = 0.0;
    double Bt = 0.0;
    double Bp = 0.0;
    for (int n_idx = 1; n_idx <= WMM_NMAX; ++n_idx)
    {
        double ar = pow(sr, n_idx + 2);
        for (int m = 0; m <= n_idx; ++m)
        {
            double g = g_coeff[n_idx][m] / 100.0;
            double h = h_coeff[n_idx][m] / 100.0;
            double cos_m = cos(m * lon_rad);
            double sin_m = sin(m * lon_rad);
            double t = g * cos_m + h * sin_m;
            Br += ar * (n_idx + 1) * t * P[n_idx][m];
            Bt -= ar * t * dP[n_idx][m];
            if (m != 0)
            {
                Bp += ar * m * (g * sin_m - h * cos_m) * P[n_idx][m] / ct;
            }
        }
    }

    // Convert spherical to NED
    *n = -Bt;
    *e = Bp;
    *d = -Br;
    return 0;
}
