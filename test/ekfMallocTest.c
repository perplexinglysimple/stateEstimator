#include "EKF.h"

// System dynamics model function prototype
void TransitionFunction(EKFMatrix* x, EKFMatrix* x_predicted, EKFState* ekf, void* userData);
// Measurement function prototype
void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);

int main()
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    // This should fail because the options struct is not initialized.
    if (EKFInit(&ekf, &options) == EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to catch uninitialized options struct.");
        return -1;
    }
    // Initialize the options struct.
    options.x0 = NULL;
    options.P0 = NULL;
    options.Q = NULL;
    options.R = NULL;
    options.A = NULL;
    options.n = 0;
    options.f = NULL;
    options.h = NULL;
    options.numberOfStates = 0;
    options.numberOfMeasurements = 0;
    options.useFiniteDifferenceJacobian = false;
    options.mallocFlag = true;
    // This should fail because the options struct has .
    if (EKFInit(&ekf, &options) == EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to catch uninitialized options struct.");
        return -1;
    }
    // Initialize the options struct.
    STATIC_MATRIX_DIRECTIVE(options.x0, 1, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 1, 1, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 1, 1, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 1, 1, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 1, 1, A);
    options.n = 1;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.numberOfStates = 1;
    options.numberOfMeasurements = 1;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;
    PRE_INIT_ALLOC(&ekf, 1, 1, true);
    // This should succeed because the options struct has been initialized.
    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to initialize the EKF.");
        return -1;
    }
    // Cleanup the EKF.
    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed to cleanup the EKF.");
        return -1;
    }
    return 0;
}

void TransitionFunction(EKFMatrix *x, EKFMatrix *x_predicted, EKFState* ekf, void* userData)
{
    (void)x;
    (void)x_predicted;
    (void)ekf;
    (void)userData;
    LOG_INFO("TransitionFunction() called.");
}

void MeasurementFunction(EKFMatrix *x, EKFMatrix *z, EKFState* ekf, void* userData)
{
    (void)x;
    (void)z;
    (void)ekf;
    (void)userData;
    LOG_INFO("MeasurementFunction() called.");
}
