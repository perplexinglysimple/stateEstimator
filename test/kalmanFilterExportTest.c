/**
 * @file kalmanFilterExportTest.c
 * @brief Export linear Kalman filter outputs for Python comparison.
 */

#include "kalmanFilter.h"
#include "matrixMath.h"
#include "utils.h"
#include <stdio.h>

int main(int argc, char** argv)
{
    const char* outPath = "kalman_filterpy_compare.csv";
    if (argc > 1 && argv[1] != NULL)
    {
        outPath = argv[1];
    }

    struct matrix* x = NULL;
    struct matrix* P = NULL;
    struct matrix* F = NULL;
    struct matrix* H = NULL;
    struct matrix* R = NULL;
    struct matrix* Q = NULL;
    struct matrix* z = NULL;

    struct matrix* F_T = NULL;
    struct matrix* H_T = NULL;
    struct matrix* temp_x = NULL;
    struct matrix* temp_fp = NULL;
    struct matrix* temp_fpft = NULL;
    struct matrix* temp_pht = NULL;

    INIT_MATRIX(x, 2, 1);
    INIT_MATRIX(P, 2, 2);
    INIT_MATRIX(F, 2, 2);
    INIT_MATRIX(H, 2, 2);
    INIT_MATRIX(R, 2, 2);
    INIT_MATRIX(Q, 2, 2);
    INIT_MATRIX(z, 2, 1);

    INIT_MATRIX(F_T, 2, 2);
    INIT_MATRIX(H_T, 2, 2);
    INIT_MATRIX(temp_x, 2, 1);
    INIT_MATRIX(temp_fp, 2, 2);
    INIT_MATRIX(temp_fpft, 2, 2);
    INIT_MATRIX(temp_pht, 2, 2);

    // Initial state and model config (match scripts/compare_kalman.py)
    SET_MATRIX(*x, 0, 0, 0.0);
    SET_MATRIX(*x, 1, 0, 1.0);

    SET_MATRIX(*P, 0, 0, 1.0);
    SET_MATRIX(*P, 0, 1, 0.0);
    SET_MATRIX(*P, 1, 0, 0.0);
    SET_MATRIX(*P, 1, 1, 1.0);

    SET_MATRIX(*F, 0, 0, 1.0);
    SET_MATRIX(*F, 0, 1, 0.1);
    SET_MATRIX(*F, 1, 0, 0.0);
    SET_MATRIX(*F, 1, 1, 1.0);

    SET_MATRIX(*H, 0, 0, 1.0);
    SET_MATRIX(*H, 0, 1, 0.0);
    SET_MATRIX(*H, 1, 0, 0.0);
    SET_MATRIX(*H, 1, 1, 1.0);

    SET_MATRIX(*Q, 0, 0, 0.01);
    SET_MATRIX(*Q, 0, 1, 0.0);
    SET_MATRIX(*Q, 1, 0, 0.0);
    SET_MATRIX(*Q, 1, 1, 0.01);

    SET_MATRIX(*R, 0, 0, 0.04);
    SET_MATRIX(*R, 0, 1, 0.0);
    SET_MATRIX(*R, 1, 0, 0.0);
    SET_MATRIX(*R, 1, 1, 0.04);

    struct kalmanMatrixes k = {
        .x_ = x,
        .P_ = P,
        .F_ = F,
        .H_ = H,
        .R_ = R,
        .Q_ = Q,
        .z_ = z,
        .F_TRANSPOSE = F_T,
        .H_TRANSPOSE = H_T,
        .TEMP_X_ = temp_x,
        .TEMP_FP_ = temp_fp,
        .TEMP_FPFT_ = temp_fpft,
        .TEMP_PHT_ = temp_pht,
    };

    if (kalmanFilterInit(&k) != MATRIX_SUCCESS)
    {
        LOG_ERROR("kalmanFilterInit() failed in export test.");
        return -1;
    }

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

    double measurements[10][2] = {
        {0.0, 1.0}, {0.1, 1.0}, {0.2, 1.0}, {0.3, 1.0}, {0.4, 1.0},
        {0.5, 1.0}, {0.6, 1.0}, {0.7, 1.0}, {0.8, 1.0}, {0.9, 1.0},
    };

    for (int step = 0; step < 10; ++step)
    {
        if (kalmanFilterPredict(&k) != MATRIX_SUCCESS)
        {
            LOG_ERROR("kalmanFilterPredict() failed in export test.");
            fclose(f);
            return -1;
        }

        SET_MATRIX(*z, 0, 0, measurements[step][0]);
        SET_MATRIX(*z, 1, 0, measurements[step][1]);

        if (kalmanFilterUpdate(&k) != MATRIX_SUCCESS)
        {
            LOG_ERROR("kalmanFilterUpdate() failed in export test.");
            fclose(f);
            return -1;
        }

        fprintf(f, "%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n", step, ACCESS_MATRIX(*x, 0, 0), ACCESS_MATRIX(*x, 1, 0),
                ACCESS_MATRIX(*P, 0, 0), ACCESS_MATRIX(*P, 0, 1), ACCESS_MATRIX(*P, 1, 0), ACCESS_MATRIX(*P, 1, 1));
    }

    fclose(f);
    FREE_MATRIX(x);
    FREE_MATRIX(P);
    FREE_MATRIX(F);
    FREE_MATRIX(H);
    FREE_MATRIX(R);
    FREE_MATRIX(Q);
    FREE_MATRIX(z);
    FREE_MATRIX(F_T);
    FREE_MATRIX(H_T);
    FREE_MATRIX(temp_x);
    FREE_MATRIX(temp_fp);
    FREE_MATRIX(temp_fpft);
    FREE_MATRIX(temp_pht);
    return 0;
}
