#include "EKF.h"
#include <math.h>

static void TransitionIdentity(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) userData;
    multMatrix(ekf->A, x, x_pred);
}

static void MeasurementIdentity(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    if (z->jaggedAlloc)
    {
        z->mat[0][0] = ACCESS_MATRIX(*x, 0, 0);
    }
    else
    {
        ACCESS_STATIC_MATRIX(*z, 0, 0) = ACCESS_MATRIX(*x, 0, 0);
    }
}

int main(void)
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};

    PRE_INIT_ALLOC(&ekf, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options.x0, 1, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 1, 1, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 1, 1, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 1, 1, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 1, 1, A);

    ACCESS_STATIC_MATRIX(*options.x0, 0, 0) = 2.0;
    ACCESS_STATIC_MATRIX(*options.P0, 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*options.Q, 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*options.R, 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*options.A, 0, 0) = 1.0;

    options.n = 1;
    options.numberOfStates = 1;
    options.numberOfMeasurements = 1;
    options.f = TransitionIdentity;
    options.h = MeasurementIdentity;
    options.updateAMatrix = NULL;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit failed.");
        return -1;
    }

    if (EKFPredict(&ekf, 0.0, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict failed with NULL updateAMatrix.");
        return -1;
    }

    if (fabs(ACCESS_MATRIX(*ekf.x, 0, 0) - 2.0) > 1e-9)
    {
        LOG_ERROR("State prediction changed unexpectedly.");
        return -1;
    }

    if (fabs(ACCESS_MATRIX(*ekf.P, 0, 0) - 1.1) > 1e-4)
    {
        LOG_ERROR("Covariance prediction mismatch.");
        return -1;
    }

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup failed.");
        return -1;
    }

    return 0;
}
