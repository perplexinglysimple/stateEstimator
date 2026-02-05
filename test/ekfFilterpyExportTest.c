#include "EKF.h"
#include <stdio.h>

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

int main(int argc, char** argv)
{
    const char* outPath = "ekf_filterpy_compare.csv";
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
        LOG_ERROR("EKFInit() failed in filterpy export test.");
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, 2, 1, z);

    // Deterministic measurements (z = [pos, vel])
    ekfType measurements[10][2] = {{0.0, 1.0}, {0.1, 1.0}, {0.2, 1.0}, {0.3, 1.0}, {0.4, 1.0},
                                   {0.5, 1.0}, {0.6, 1.0}, {0.7, 1.0}, {0.8, 1.0}, {0.9, 1.0}};

    FILE* f = NULL;
#ifdef _WIN32
    if (fopen_s(&f, outPath, "w") != 0)
    {
        f = NULL;
    }
#else
    f = fopen(outPath, "w");
#endif
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
            LOG_ERROR("EKFPredict() failed in export test.");
            fclose(f);
            return -1;
        }
        ACCESS_STATIC_MATRIX(*(measurement.z), 0, 0) = measurements[k][0];
        ACCESS_STATIC_MATRIX(*(measurement.z), 1, 0) = measurements[k][1];
        if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFUpdate() failed in export test.");
            fclose(f);
            return -1;
        }

        fprintf(f, "%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", k, ACCESS_MATRIX(*(ekf.x), 0, 0),
                ACCESS_MATRIX(*(ekf.x), 1, 0), ACCESS_MATRIX(*(ekf.P), 0, 0), ACCESS_MATRIX(*(ekf.P), 0, 1),
                ACCESS_MATRIX(*(ekf.P), 1, 0), ACCESS_MATRIX(*(ekf.P), 1, 1));
    }

    fclose(f);
    EKFCleanup(&ekf);
    return 0;
}

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) userData;
    multMatrix(ekf->A, x, x_pred);
}

static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    copyMatrix(x, z);
}

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void) x;
    (void) userData;
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
    (void) x;
    (void) ekf;
    (void) userData;
    for (int i = 0; i < J->row; ++i)
    {
        for (int j = 0; j < J->col; ++j)
        {
            ekfType val = (i == j) ? 1 : 0;
            if (J->jaggedAlloc)
            {
                J->mat[i][j] = val;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*J, i, j) = val;
            }
        }
    }
}

static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    (void) x;
    (void) ekf;
    (void) time;
    (void) userData;
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
