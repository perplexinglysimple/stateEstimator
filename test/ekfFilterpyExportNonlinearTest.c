#include "EKF.h"
#include <stdio.h>
#include <math.h>

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

int main(int argc, char** argv)
{
    const char* outPath = "ekf_filterpy_compare_nonlinear.csv";
    if (argc > 1 && argv[1] != NULL)
    {
        outPath = argv[1];
    }

    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    PRE_INIT_ALLOC(&ekf, 2, 2, true);

    STATIC_MATRIX_DIRECTIVE(options.x0, 2, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 2, 2, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 2, 2, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 2, 2, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 2, 2, A);

    // Model config
    ACCESS_STATIC_MATRIX(*(options.x0), 0, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.x0), 1, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.P0), 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.P0), 1, 1) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.Q), 0, 0) = 0.01;
    ACCESS_STATIC_MATRIX(*(options.Q), 1, 1) = 0.01;
    ACCESS_STATIC_MATRIX(*(options.R), 0, 0) = 0.04;
    ACCESS_STATIC_MATRIX(*(options.R), 1, 1) = 0.04;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 0) = 1.0;
    ACCESS_STATIC_MATRIX(*(options.A), 0, 1) = 0.1;
    ACCESS_STATIC_MATRIX(*(options.A), 1, 0) = 0.0;
    ACCESS_STATIC_MATRIX(*(options.A), 1, 1) = 1.0;

    options.n = 2;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.updateAMatrix = UpdateAMatrix;
    options.numberOfStates = 2;
    options.numberOfMeasurements = 2;
    options.useFiniteDifferenceJacobian = false;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in filterpy nonlinear export test.");
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 2, 1, z);

    // Deterministic measurements (z = [pos, vel]) for nonlinear compare
    ekfType measurements[10][2] = {
        {0.00, 1.00},
        {0.05, 0.98},
        {0.10, 0.97},
        {0.16, 0.95},
        {0.21, 0.94},
        {0.27, 0.92},
        {0.33, 0.90},
        {0.39, 0.89},
        {0.46, 0.87},
        {0.52, 0.86}
    };

    FILE* f = fopen(outPath, "w");
    if (!f)
    {
        LOG_ERROR("Failed to open output file for writing.");
        return -1;
    }
    fprintf(f, "step,x0,x1,P00,P01,P10,P11\n");

    for (int k = 0; k < 10; ++k)
    {
        if (EKFPredict(&ekf, 0.1, NULL) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFPredict() failed in nonlinear export test.");
            fclose(f);
            return -1;
        }
        ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = measurements[k][0];
        ACCESS_STATIC_MATRIX(*(measurement.z), 1, 0) = measurements[k][1];
        if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFUpdate() failed in nonlinear export test.");
            fclose(f);
            return -1;
        }

        fprintf(f, "%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                k,
                ACCESS_MATRIX(*(ekf.x), 0, 0),
                ACCESS_MATRIX(*(ekf.x), 1, 0),
                ACCESS_MATRIX(*(ekf.P), 0, 0),
                ACCESS_MATRIX(*(ekf.P), 0, 1),
                ACCESS_MATRIX(*(ekf.P), 1, 0),
                ACCESS_MATRIX(*(ekf.P), 1, 1));
    }

    fclose(f);
    EKFCleanup(&ekf);
    return 0;
}

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    double dt = 0.1;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double y0 = x0 + x1 * dt + 0.5 * dt * dt * sin(x0);
    double y1 = x1 + dt * cos(x1);
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

static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
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

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    double dt = 0.1;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double f00 = 1.0 + 0.5 * dt * dt * cos(x0);
    double f01 = dt;
    double f10 = 0.0;
    double f11 = 1.0 - dt * sin(x1);

    if (J->jaggedAlloc)
    {
        J->mat[0][0] = f00;
        J->mat[0][1] = f01;
        J->mat[1][0] = f10;
        J->mat[1][1] = f11;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*J, 0, 0) = f00;
        ACCESS_STATIC_MATRIX(*J, 0, 1) = f01;
        ACCESS_STATIC_MATRIX(*J, 1, 0) = f10;
        ACCESS_STATIC_MATRIX(*J, 1, 1) = f11;
    }
}

static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    double x0 = ACCESS_MATRIX(*x, 0, 0);
    double x1 = ACCESS_MATRIX(*x, 1, 0);
    double h00 = 2.0 * x0;
    double h01 = 0.0;
    double h10 = 0.0;
    double h11 = 2.0 * x1;

    if (J->jaggedAlloc)
    {
        J->mat[0][0] = h00;
        J->mat[0][1] = h01;
        J->mat[1][0] = h10;
        J->mat[1][1] = h11;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*J, 0, 0) = h00;
        ACCESS_STATIC_MATRIX(*J, 0, 1) = h01;
        ACCESS_STATIC_MATRIX(*J, 1, 0) = h10;
        ACCESS_STATIC_MATRIX(*J, 1, 1) = h11;
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
        A->mat[0][1] = 0.1;
        A->mat[1][0] = 0.0;
        A->mat[1][1] = 1.0;
    }
    else
    {
        ACCESS_STATIC_MATRIX(*A, 0, 0) = 1.0;
        ACCESS_STATIC_MATRIX(*A, 0, 1) = 0.1;
        ACCESS_STATIC_MATRIX(*A, 1, 0) = 0.0;
        ACCESS_STATIC_MATRIX(*A, 1, 1) = 1.0;
    }
}
