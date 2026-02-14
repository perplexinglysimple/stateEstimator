#include "EKF.h"
#include <math.h>

static void TransitionFD(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void MeasurementFD(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void IdentityUpdateA(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

static int assertNear(double actual, double expected, double eps, const char* msg)
{
    (void) msg;
    if (fabs(actual - expected) > eps)
    {
        LOG_ERROR("%s: expected %.6f got %.6f", msg, expected, actual);
        return -1;
    }
    return 0;
}

int main()
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    PRE_INIT_ALLOC(&ekf, 2, 2, true);

    STATIC_MATRIX_DIRECTIVE(options.x0, 2, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 2, 2, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 2, 2, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 2, 2, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 2, 2, A);

    ACCESS_STATIC_MATRIX(*(options.x0), 0, 0) = 0.5;
    ACCESS_STATIC_MATRIX(*(options.x0), 1, 0) = -0.4;
    ACCESS_STATIC_MATRIX(*(options.P0), 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*(options.P0), 1, 1) = 1;
    ACCESS_STATIC_MATRIX(*(options.Q), 0, 0) = 0;
    ACCESS_STATIC_MATRIX(*(options.Q), 1, 1) = 0;
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 0.01;
    ACCESS_STATIC_MATRIX(*(options.R), 1, 1) = 0.01;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*(options.A), 1, 1) = 1;

    options.n = 2;
    options.f = TransitionFD;
    options.h = MeasurementFD;
    options.updateAMatrix = IdentityUpdateA;
    options.numberOfStates = 2;
    options.numberOfMeasurements = 2;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in finite-diff test.");
        return -1;
    }

    if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict() failed in finite-diff test.");
        return -1;
    }

    // Expected Jacobian for f(x) evaluated at predicted state:
    // f0 = x0^2 + x1
    // f1 = sin(x0) + cos(x1)
    // F = [[2*x0, 1], [cos(x0), -sin(x1)]]
    double x0 = ACCESS_MATRIX(*(ekf.x), 0, 0);
    double x1 = ACCESS_MATRIX(*(ekf.x), 1, 0);
    double f00 = 2.0 * x0;
    double f01 = 1.0;
    double f10 = cos(x0);
    double f11 = -sin(x1);

    if (assertNear(ACCESS_MATRIX(*(ekf._F), 0, 0), f00, 1e-3, "F(0,0)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._F), 0, 1), f01, 1e-3, "F(0,1)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._F), 1, 0), f10, 1e-3, "F(1,0)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._F), 1, 1), f11, 1e-3, "F(1,1)") != 0)
        return -1;

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 2, 1, z);
    ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*(measurement.z), 1, 0) = -0.2;

    if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdate() failed in finite-diff test.");
        return -1;
    }

    // Measurement h(x) evaluated at predicted state (_x_predicted):
    // h0 = x0^2
    // h1 = x1^2
    // H = [[2*x0, 0], [0, 2*x1]]
    double x0h = ACCESS_MATRIX(*(ekf._x_predicted), 0, 0);
    double x1h = ACCESS_MATRIX(*(ekf._x_predicted), 1, 0);
    double h00 = 2.0 * x0h;
    double h01 = 0.0;
    double h10 = 0.0;
    double h11 = 2.0 * x1h;

    if (assertNear(ACCESS_MATRIX(*(ekf._H), 0, 0), h00, 1e-3, "H(0,0)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._H), 0, 1), h01, 1e-3, "H(0,1)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._H), 1, 0), h10, 1e-3, "H(1,0)") != 0)
        return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf._H), 1, 1), h11, 1e-3, "H(1,1)") != 0)
        return -1;

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed in finite-diff test.");
        return -1;
    }

    return 0;
}

static void TransitionFD(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double y0 = x0 * x0 + x1;
    double y1 = sin(x0) + cos(x1);
    if (x_pred->jaggedAlloc)
    {
        x_pred->mat[0][0] = y0;
        x_pred->mat[1][0] = y1;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*x_pred, 0, 0) = y0;
        ACCESS_STATIC_MATRIX(*x_pred, 1, 0) = y1;
    }
}

static void MeasurementFD(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double z0 = x0 * x0;
    double z1 = x1 * x1;
    if (z->jaggedAlloc)
    {
        z->mat[0][0] = z0;
        z->mat[1][0] = z1;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*z, 0, 0) = z0;
        ACCESS_STATIC_MATRIX(*z, 1, 0) = z1;
    }
}

static void IdentityUpdateA(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    (void) x;
    (void) ekf;
    (void) time;
    (void) userData;
    if (A->jaggedAlloc)
    {
        A->mat[0][0] = 1;
        A->mat[0][1] = 0;
        A->mat[1][0] = 0;
        A->mat[1][1] = 1;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*A, 0, 0) = 1;
        ACCESS_STATIC_MATRIX(*A, 0, 1) = 0;
        ACCESS_STATIC_MATRIX(*A, 1, 0) = 0;
        ACCESS_STATIC_MATRIX(*A, 1, 1) = 1;
    }
}
