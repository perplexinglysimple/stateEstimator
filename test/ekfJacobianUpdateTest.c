#include "EKF.h"
#include <math.h>

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_predicted, EKFState* ekf, void* userData);
static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
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
    // 1) Missing Jacobian callbacks should fail when finite diff is disabled
    EKFState ekf_missing = {0};
    EKFConfigOptions options_missing = {0};
    PRE_INIT_ALLOC(&ekf_missing, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options_missing.x0, 1, 1, x0m);
    STATIC_MATRIX_DIRECTIVE(options_missing.P0, 1, 1, P0m);
    STATIC_MATRIX_DIRECTIVE(options_missing.Q, 1, 1, Qm);
    STATIC_MATRIX_DIRECTIVE(options_missing.R, 1, 1, Rm);
    STATIC_MATRIX_DIRECTIVE(options_missing.A, 1, 1, Am);
    options_missing.n = 1;
    options_missing.f = TransitionFunction;
    options_missing.h = MeasurementFunction;
    options_missing.updateAMatrix = UpdateAMatrix;
    options_missing.numberOfStates = 1;
    options_missing.numberOfMeasurements = 1;
    options_missing.mallocFlag = true;
    options_missing.useFiniteDifferenceJacobian = false;
    options_missing.jacobianF = NULL;
    options_missing.jacobianH = NULL;
    if (EKFInit(&ekf_missing, &options_missing) == EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() should fail when Jacobian callbacks are missing.");
        return -1;
    }
    EKFCleanup(&ekf_missing);

    // 2) Callbacks provided should succeed
    EKFState ekf_callbacks = {0};
    EKFConfigOptions options_callbacks = {0};
    PRE_INIT_ALLOC(&ekf_callbacks, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options_callbacks.x0, 1, 1, x0c);
    STATIC_MATRIX_DIRECTIVE(options_callbacks.P0, 1, 1, P0c);
    STATIC_MATRIX_DIRECTIVE(options_callbacks.Q, 1, 1, Qc);
    STATIC_MATRIX_DIRECTIVE(options_callbacks.R, 1, 1, Rc);
    STATIC_MATRIX_DIRECTIVE(options_callbacks.A, 1, 1, Ac);
    options_callbacks.n = 1;
    options_callbacks.f = TransitionFunction;
    options_callbacks.h = MeasurementFunction;
    options_callbacks.updateAMatrix = UpdateAMatrix;
    options_callbacks.numberOfStates = 1;
    options_callbacks.numberOfMeasurements = 1;
    options_callbacks.mallocFlag = true;
    options_callbacks.useFiniteDifferenceJacobian = false;
    options_callbacks.jacobianF = StateJacobianFunction;
    options_callbacks.jacobianH = MeasurementJacobianFunction;
    if (EKFInit(&ekf_callbacks, &options_callbacks) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed with Jacobian callbacks provided.");
        return -1;
    }
    EKFCleanup(&ekf_callbacks);

    // 3) Finite difference should succeed without callbacks
    EKFState ekf_fd = {0};
    EKFConfigOptions options_fd = {0};
    PRE_INIT_ALLOC(&ekf_fd, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options_fd.x0, 1, 1, x0f);
    STATIC_MATRIX_DIRECTIVE(options_fd.P0, 1, 1, P0f);
    STATIC_MATRIX_DIRECTIVE(options_fd.Q, 1, 1, Qf);
    STATIC_MATRIX_DIRECTIVE(options_fd.R, 1, 1, Rf);
    STATIC_MATRIX_DIRECTIVE(options_fd.A, 1, 1, Af);
    options_fd.n = 1;
    options_fd.f = TransitionFunction;
    options_fd.h = MeasurementFunction;
    options_fd.updateAMatrix = UpdateAMatrix;
    options_fd.numberOfStates = 1;
    options_fd.numberOfMeasurements = 1;
    options_fd.mallocFlag = true;
    options_fd.useFiniteDifferenceJacobian = true;
    options_fd.jacobianF = NULL;
    options_fd.jacobianH = NULL;
    if (EKFInit(&ekf_fd, &options_fd) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed with finite difference Jacobian enabled.");
        return -1;
    }
    EKFCleanup(&ekf_fd);

    // 4) Predict/update correctness on 1D linear model
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    PRE_INIT_ALLOC(&ekf, 1, 1, true);
    STATIC_MATRIX_DIRECTIVE(options.x0, 1, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 1, 1, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 1, 1, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 1, 1, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 1, 1, A);
    ACCESS_STATIC_MATRIX(*(options.x0), 0, 0) = 0;
    ACCESS_STATIC_MATRIX(*(options.P0), 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*(options.Q), 0, 0) = 0;
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 0) = 1;
    options.n = 1;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.updateAMatrix = UpdateAMatrix;
    options.numberOfStates = 1;
    options.numberOfMeasurements = 1;
    options.mallocFlag = true;
    options.useFiniteDifferenceJacobian = false;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to initialize EKF for correctness test.");
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 1, 1, z);
    ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = 1;

    if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFPredict() failed in correctness test.");
        return -1;
    }

    double p_before = ACCESS_MATRIX(*(ekf.P), 0, 0);
    double z_meas = ACCESS_MATRIX(*(measurement.z), 0, 0);
    double x_pred = ACCESS_MATRIX(*(ekf._x_predicted), 0, 0);
    if (x_pred != 0)
    {
        LOG_ERROR("Predicted state expected 0 but got %.6f", x_pred);
        return -1;
    }
    if (z_meas != 1)
    {
        LOG_ERROR("Measurement z was not set correctly.");
        return -1;
    }

    if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFUpdate() failed in correctness test.");
        return -1;
    }

    double x_est = ACCESS_MATRIX(*(ekf.x), 0, 0);
    double p_after = ACCESS_MATRIX(*(ekf.P), 0, 0);
    double k_val = ACCESS_MATRIX(*(ekf._K), 0, 0);
    if (k_val == 0)
    {
        LOG_ERROR("Kalman gain K is zero in correctness test.");
        return -1;
    }

    if (assertNear(x_est, 0.9090909, 1e-4, "x estimate") != 0)
    {
        return -1;
    }
    if (assertNear(p_after, 0.0909091, 1e-4, "P update") != 0)
    {
        return -1;
    }
    if (p_after >= p_before)
    {
        LOG_ERROR("P should decrease after measurement update.");
        return -1;
    }

    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed.");
        return -1;
    }

    return 0;
}

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_predicted, EKFState* ekf, void* userData)
{
    (void)userData;
    multMatrix(ekf->A, x, x_predicted);
}

static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    copyMatrix(x, z);
}

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)userData;
    if (J->jaggedAlloc)
    {
        J->mat[0][0] = ACCESS_MATRIX(*(ekf->A), 0, 0);
    }
    else
    {
        ACCESS_STATIC_MATRIX(*J, 0, 0) = ACCESS_MATRIX(*(ekf->A), 0, 0);
    }
}

static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)ekf;
    (void)userData;
    if (J->jaggedAlloc)
    {
        J->mat[0][0] = 1;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*J, 0, 0) = 1;
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
        A->mat[0][0] = 1;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*A, 0, 0) = 1;
    }
}
