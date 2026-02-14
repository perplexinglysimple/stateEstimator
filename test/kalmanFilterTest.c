/**
 * @file kalmanFilterTest.c
 */

#include "kalmanFilter.h"
#include "matrixMath.h"
#include "utils.h"
#include <math.h>

// Tests return matrixReturnCodes; override EKF macro to avoid enum-conversion warnings.
#ifdef MATRIX_MATH_RETURN_CHECK
#undef MATRIX_MATH_RETURN_CHECK
#endif
#define MATRIX_MATH_RETURN_CHECK(ret)                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        matrixReturnCodes _ret = (ret);                                                                                \
        if (_ret != MATRIX_SUCCESS)                                                                                    \
        {                                                                                                              \
            return _ret;                                                                                               \
        }                                                                                                              \
    } while (0)

static matrixReturnCodes expectReturnCode(matrixReturnCodes got, matrixReturnCodes expected, const char* msg)
{
    (void) msg;
    if (got != expected)
    {
        LOG_ERROR("%s (got %d, expected %d)", msg, got, expected);
        return MATRIX_ERROR;
    }
    return MATRIX_SUCCESS;
}

static matrixReturnCodes expectScalarApprox(matrixType got, matrixType expected, matrixType eps, const char* msg)
{
    (void) msg;
    if (fabs((double) (got - expected)) > (double) eps)
    {
        LOG_ERROR("%s (got %.10f, expected %.10f)", msg, (double) got, (double) expected);
        return MATRIX_ERROR;
    }
    return MATRIX_SUCCESS;
}

static matrixReturnCodes kalmanLinear1DTest(void)
{
    LOG_INFO("Starting Kalman 1D Test");

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

    INIT_MATRIX(x, 1, 1);
    INIT_MATRIX(P, 1, 1);
    INIT_MATRIX(F, 1, 1);
    INIT_MATRIX(H, 1, 1);
    INIT_MATRIX(R, 1, 1);
    INIT_MATRIX(Q, 1, 1);
    INIT_MATRIX(z, 1, 1);

    INIT_MATRIX(F_T, 1, 1);
    INIT_MATRIX(H_T, 1, 1);
    INIT_MATRIX(temp_x, 1, 1);
    INIT_MATRIX(temp_fp, 1, 1);
    INIT_MATRIX(temp_fpft, 1, 1);
    INIT_MATRIX(temp_pht, 1, 1);

    SET_MATRIX(*x, 0, 0, 0);
    SET_MATRIX(*P, 0, 0, 1);
    SET_MATRIX(*F, 0, 0, 1);
    SET_MATRIX(*H, 0, 0, 1);
    SET_MATRIX(*R, 0, 0, (matrixType) 0.2);
    SET_MATRIX(*Q, 0, 0, (matrixType) 0.1);
    SET_MATRIX(*z, 0, 0, 1);

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

    MATRIX_MATH_RETURN_CHECK(kalmanFilterInit(&k));

    MATRIX_MATH_RETURN_CHECK(kalmanFilterPredict(&k));
    MATRIX_MATH_RETURN_CHECK(expectScalarApprox(ACCESS_MATRIX(*x, 0, 0), 0.0, (matrixType) 1e-9, "predict x"));
    MATRIX_MATH_RETURN_CHECK(expectScalarApprox(ACCESS_MATRIX(*P, 0, 0), 1.1, (matrixType) 1e-9, "predict P"));

    MATRIX_MATH_RETURN_CHECK(kalmanFilterUpdate(&k));
    MATRIX_MATH_RETURN_CHECK(expectScalarApprox(ACCESS_MATRIX(*x, 0, 0), 0.846153846, (matrixType) 1e-6, "update x"));
    MATRIX_MATH_RETURN_CHECK(expectScalarApprox(ACCESS_MATRIX(*P, 0, 0), 0.169230769, (matrixType) 1e-6, "update P"));

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

    LOG_INFO("Completed Kalman 1D Test");
    return MATRIX_SUCCESS;
}

static matrixReturnCodes kalmanInitDimensionMismatchTest(void)
{
    LOG_INFO("Starting Kalman Init Dimension Mismatch Test");

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
    INIT_MATRIX(H, 1, 2); // mismatch: measurement dim != state dim
    INIT_MATRIX(R, 1, 1);
    INIT_MATRIX(Q, 2, 2);
    INIT_MATRIX(z, 1, 1);

    INIT_MATRIX(F_T, 2, 2);
    INIT_MATRIX(H_T, 2, 1);
    INIT_MATRIX(temp_x, 2, 1);
    INIT_MATRIX(temp_fp, 2, 2);
    INIT_MATRIX(temp_fpft, 2, 2);
    INIT_MATRIX(temp_pht, 2, 1);

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

    MATRIX_MATH_RETURN_CHECK(expectReturnCode(kalmanFilterInit(&k), MATRIX_DIMENSION_MISMATCH, "init mismatch"));

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

    LOG_INFO("Completed Kalman Init Dimension Mismatch Test");
    return MATRIX_SUCCESS;
}

int main(void)
{
    LOG_INFO("Starting Kalman Filter Tests");

    MATRIX_MATH_RETURN_CHECK(kalmanLinear1DTest());
    MATRIX_MATH_RETURN_CHECK(kalmanInitDimensionMismatchTest());

    LOG_INFO("Completed Kalman Filter Tests");
    return 0;
}
