#include "EKF.h"
#include <math.h>

static void TransitionIdentity(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void MeasurementX(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

static int assertNear(double actual, double expected, double eps, const char* msg)
{
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
    PRE_INIT_ALLOC(&ekf, 2, 1, true);

    STATIC_MATRIX_DIRECTIVE(options.x0, 2, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 2, 2, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 2, 2, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 1, 1, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 2, 2, A);

    ACCESS_STATIC_MATRIX(*(options.x0), 0, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.x0), 1, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.P0), 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.P0), 1, 1) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.Q), 0, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.Q), 1, 1) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.A), 1, 1) = 1.0;

    options.n = 2;
    options.f = TransitionIdentity;
    options.h = MeasurementX;
    options.updateAMatrix = UpdateAMatrix;
    options.numberOfStates = 2;
    options.numberOfMeasurements = 1;
    options.useFiniteDifferenceJacobian = false;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in measurement dimension test.");
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 1, 1, z);
    ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = 1.0;

    if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict() failed in measurement dimension test.");
        return -1;
    }

    if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdate() failed in measurement dimension test.");
        return -1;
    }

    if (assertNear(ACCESS_MATRIX(*(ekf.x), 0, 0), 0.9090909, 1e-4, "x0 update") != 0) return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf.x), 1, 0), 0.0, 1e-6, "x1 update") != 0) return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf.P), 0, 0), 0.0909091, 1e-4, "P00 update") != 0) return -1;
    if (assertNear(ACCESS_MATRIX(*(ekf.P), 1, 1), 1.0, 1e-6, "P11 update") != 0) return -1;

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed in measurement dimension test.");
        return -1;
    }

    return 0;
}

static void TransitionIdentity(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void)userData;
    multMatrix(ekf->A, x, x_pred);
}

static void MeasurementX(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    if (z->jaggedAlloc)
    {
        z->mat[0][0] = ACCESS_MATRIX(*x, 0, 0);
    }
    else
    {
        ACCESS_STATIC_MATRIX(*z, 0, 0) = ACCESS_MATRIX(*x, 0, 0);
    }
}

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)userData;
    for (int i = 0; i < J->row; ++i)
    {
        for (int j = 0; j < J->col; ++j)
        {
            if (J->jaggedAlloc)
            {
                J->mat[i][j] = ACCESS_MATRIX(*(ekf->A), i, j);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*J, i, j) = ACCESS_MATRIX(*(ekf->A), i, j);
            }
        }
    }
}

static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)ekf;
    (void)userData;
    if (J->jaggedAlloc)
    {
        J->mat[0][0] = 1.0;
        J->mat[0][1] = 0.0;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*J, 0, 0) = 1.0;
        ACCESS_STATIC_MATRIX(*J, 0, 1) = 0.0;
    }
}

static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    (void)x;
    (void)ekf;
    (void)time;
    (void)userData;
    if (A->jaggedAlloc)
    {
        A->mat[0][0] = 1.0;
        A->mat[0][1] = 0.0;
        A->mat[1][0] = 0.0;
        A->mat[1][1] = 1.0;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*A, 0, 0) = 1.0;
        ACCESS_STATIC_MATRIX(*A, 0, 1) = 0.0;
        ACCESS_STATIC_MATRIX(*A, 1, 0) = 0.0;
        ACCESS_STATIC_MATRIX(*A, 1, 1) = 1.0;
    }
}
