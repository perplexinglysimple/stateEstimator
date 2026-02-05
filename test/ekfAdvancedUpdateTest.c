#include "EKF.h"
#include <math.h>

static void NonlinearTransition(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void NonlinearMeasurement(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void IdentityUpdateA(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

static int matrixIsSymmetric(EKFMatrix* m, double eps)
{
    if (m->row != m->col)
    {
        return 0;
    }
    for (int i = 0; i < m->row; ++i)
    {
        for (int j = i + 1; j < m->col; ++j)
        {
            double a = ACCESS_MATRIX(*m, i, j);
            double b = ACCESS_MATRIX(*m, j, i);
            if (fabs(a - b) > eps)
            {
                return 0;
            }
        }
    }
    return 1;
}

static int matrixDiagonalPositive(EKFMatrix* m, double eps)
{
    if (m->row != m->col)
    {
        return 0;
    }
    for (int i = 0; i < m->row; ++i)
    {
        if (ACCESS_MATRIX(*m, i, i) <= eps)
        {
            return 0;
        }
    }
    return 1;
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

    // Setup: simple nonlinear transition and measurement
    ACCESS_STATIC_MATRIX(*(options.x0), 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.x0), 1, 0) = -0.2;
    ACCESS_STATIC_MATRIX(*(options.P0), 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*(options.P0), 1, 1) = 1;
    ACCESS_STATIC_MATRIX(*(options.Q), 0, 0) = 0.001;
    ACCESS_STATIC_MATRIX(*(options.Q), 1, 1) = 0.001;
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.R), 1, 1) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*(options.A), 1, 1) = 1;

    options.n = 2;
    options.f = NonlinearTransition;
    options.h = NonlinearMeasurement;
    options.updateAMatrix = IdentityUpdateA;
    options.numberOfStates = 2;
    options.numberOfMeasurements = 2;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in advanced test.");
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 2, 1, z);
    ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = 0.05;
    ACCESS_STATIC_MATRIX(*(measurement.z), 1, 0) = -0.1;

    if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict() failed in advanced test.");
        return -1;
    }

    if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdate() failed in advanced test.");
        return -1;
    }

    if (!matrixIsSymmetric(ekf.P, 1e-6))
    {
        LOG_ERROR("Joseph form covariance update should remain symmetric.");
        return -1;
    }

    if (!matrixDiagonalPositive(ekf.P, 0.0))
    {
        LOG_ERROR("Covariance diagonal should remain positive after update.");
        return -1;
    }

    // Jitter handling: create near-singular S by setting R extremely small
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 1e-12;
    ACCESS_STATIC_MATRIX(*(options.R), 1, 1) = 1e-12;
    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed before jitter test.");
        return -1;
    }

    if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict() failed before jitter test.");
        return -1;
    }
    if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdate() failed during jitter test.");
        return -1;
    }

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed in advanced test.");
        return -1;
    }

    return 0;
}

static void NonlinearTransition(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double y0 = x0 + 0.1 * sin(x1);
    double y1 = x1 + 0.1 * cos(x0);
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

static void NonlinearMeasurement(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
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
