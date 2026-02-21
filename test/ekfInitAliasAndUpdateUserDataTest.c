#include "EKF.h"
#include <math.h>

typedef struct MeasurementParams_
{
    double offset;
} MeasurementParams;

static void TransitionIdentity(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) userData;
    multMatrix(ekf->A, x, x_pred);
}

static void MeasurementWithOffset(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    MeasurementParams* params = (MeasurementParams*) userData;
    double offset = 0.0;

    (void) ekf;
    if (params != NULL)
    {
        offset = params->offset;
    }

    if (z->jaggedAlloc)
    {
        z->mat[0][0] = ACCESS_MATRIX(*x, 0, 0) + offset;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*z, 0, 0) = ACCESS_MATRIX(*x, 0, 0) + offset;
    }
}

static int nearlyEqual(double a, double b, double eps)
{
    return fabs(a - b) <= eps;
}

int main(void)
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    EKFMeasurement measurement = {0};
    EKFMeasurement badMeasurement = {0};
    MeasurementParams params;

    params.offset = 1.5;

    PRE_INIT_ALLOC(&ekf, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options.x0, 1, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 1, 1, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 1, 1, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 1, 1, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 1, 1, A);
    STATIC_MATRIX_DIRECTIVE(measurement.z, 1, 1, z);

    ACCESS_STATIC_MATRIX(*options.x0, 0, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*options.P0, 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*options.Q, 0, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*options.R, 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*options.A, 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*measurement.z, 0, 0) = 2.0;

    options.n = 0;
    options.numberOfStates = 1;
    options.numberOfMeasurements = 0;
    options.f = TransitionIdentity;
    options.h = MeasurementWithOffset;
    options.updateAMatrix = NULL;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit failed to honor numberOfStates alias.");
        return -1;
    }

    if (options.n != 1 || options.numberOfStates != 1 || options.numberOfMeasurements != 1)
    {
        LOG_ERROR("EKFInit did not normalize dimensions.");
        return -1;
    }

    if (EKFPredict(&ekf, 0.0, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict failed.");
        return -1;
    }

    if (EKFUpdateWithUserData(&ekf, &badMeasurement, &params) != EKF_NULL_POINTER)
    {
        LOG_ERROR("EKFUpdateWithUserData should reject NULL measurement->z.");
        return -1;
    }

    if (EKFUpdateWithUserData(&ekf, &measurement, &params) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdateWithUserData failed.");
        return -1;
    }

    if (!nearlyEqual(ACCESS_MATRIX(*ekf.x, 0, 0), 0.25, 1e-6))
    {
        LOG_ERROR("Unexpected state update with userData-aware measurement model.");
        return -1;
    }

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup failed.");
        return -1;
    }

    return 0;
}
