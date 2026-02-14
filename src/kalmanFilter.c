#include "kalmanFilter.h"
#include "matrixMath.h"
#include "utils.h"

static inline matrixReturnCodes kalmanMatrixNullCheck(struct matrix* ptr)
{
    NULL_CHECK_MATRIX(ptr);
    return MATRIX_SUCCESS;
}

#define KALMAN_NULL_CHECK_MATRIX(ptr)                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        matrixReturnCodes _mret = kalmanMatrixNullCheck((ptr));                                                        \
        if (_mret != MATRIX_SUCCESS)                                                                                   \
        {                                                                                                              \
            return _mret;                                                                                              \
        }                                                                                                              \
    } while (0)

#define KALMAN_RETURN_CHECK(ret)                                                                                       \
    do                                                                                                                 \
    {                                                                                                                  \
        matrixReturnCodes _ret = (ret);                                                                                \
        if (_ret != MATRIX_SUCCESS)                                                                                    \
        {                                                                                                              \
            return _ret;                                                                                               \
        }                                                                                                              \
    } while (0)

static matrixReturnCodes validateKalmanDimensions(const struct kalmanMatrixes* k)
{
    int n = k->x_->row;
    if (k->x_->col != 1)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }

    // This linear Kalman implementation assumes measurement dimension matches state dimension.
    if (k->H_->row != n || k->H_->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }

    if (k->P_->row != n || k->P_->col != n || k->F_->row != n || k->F_->col != n || k->Q_->row != n || k->Q_->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->R_->row != n || k->R_->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->F_TRANSPOSE->row != n || k->F_TRANSPOSE->col != n || k->H_TRANSPOSE->row != n || k->H_TRANSPOSE->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->TEMP_X_->row != n || k->TEMP_X_->col != 1)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->TEMP_FP_->row != n || k->TEMP_FP_->col != n || k->TEMP_FPFT_->row != n || k->TEMP_FPFT_->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->TEMP_PHT_->row != n || k->TEMP_PHT_->col != n)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (k->z_->row != n || k->z_->col != 1)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }

    return MATRIX_SUCCESS;
}

int kalmanFilterInit(struct kalmanMatrixes* k)
{
    if (k == NULL)
    {
        return MATRIX_NULL_POINTER;
    }

    KALMAN_NULL_CHECK_MATRIX(k->x_);
    KALMAN_NULL_CHECK_MATRIX(k->P_);
    KALMAN_NULL_CHECK_MATRIX(k->F_);
    KALMAN_NULL_CHECK_MATRIX(k->H_);
    KALMAN_NULL_CHECK_MATRIX(k->R_);
    KALMAN_NULL_CHECK_MATRIX(k->Q_);
    KALMAN_NULL_CHECK_MATRIX(k->z_);

    KALMAN_NULL_CHECK_MATRIX(k->F_TRANSPOSE);
    KALMAN_NULL_CHECK_MATRIX(k->H_TRANSPOSE);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_X_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FP_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FPFT_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_PHT_);

    NON_INIT_CHECK_MATRIX(k->x_);
    NON_INIT_CHECK_MATRIX(k->P_);
    NON_INIT_CHECK_MATRIX(k->F_);
    NON_INIT_CHECK_MATRIX(k->H_);
    NON_INIT_CHECK_MATRIX(k->R_);
    NON_INIT_CHECK_MATRIX(k->Q_);
    NON_INIT_CHECK_MATRIX(k->z_);

    NON_INIT_CHECK_MATRIX(k->F_TRANSPOSE);
    NON_INIT_CHECK_MATRIX(k->H_TRANSPOSE);
    NON_INIT_CHECK_MATRIX(k->TEMP_X_);
    NON_INIT_CHECK_MATRIX(k->TEMP_FP_);
    NON_INIT_CHECK_MATRIX(k->TEMP_FPFT_);
    NON_INIT_CHECK_MATRIX(k->TEMP_PHT_);

    return validateKalmanDimensions(k);
}

int kalmanFilterPredict(struct kalmanMatrixes* k)
{
    if (k == NULL)
    {
        return MATRIX_NULL_POINTER;
    }

    KALMAN_NULL_CHECK_MATRIX(k->x_);
    KALMAN_NULL_CHECK_MATRIX(k->P_);
    KALMAN_NULL_CHECK_MATRIX(k->F_);
    KALMAN_NULL_CHECK_MATRIX(k->Q_);
    KALMAN_NULL_CHECK_MATRIX(k->F_TRANSPOSE);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_X_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FP_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FPFT_);

    KALMAN_RETURN_CHECK(multMatrix(k->F_, k->x_, k->TEMP_X_));
    KALMAN_RETURN_CHECK(copyMatrix(k->TEMP_X_, k->x_));

    KALMAN_RETURN_CHECK(multMatrix(k->F_, k->P_, k->TEMP_FP_));
    KALMAN_RETURN_CHECK(transposeMatrix(k->F_, k->F_TRANSPOSE));
    KALMAN_RETURN_CHECK(multMatrix(k->TEMP_FP_, k->F_TRANSPOSE, k->TEMP_FPFT_));
    KALMAN_RETURN_CHECK(addMatrix(k->TEMP_FPFT_, k->Q_, k->TEMP_FP_));
    KALMAN_RETURN_CHECK(copyMatrix(k->TEMP_FP_, k->P_));

    return MATRIX_SUCCESS;
}

int kalmanFilterUpdate(struct kalmanMatrixes* k)
{
    if (k == NULL)
    {
        return MATRIX_NULL_POINTER;
    }

    KALMAN_NULL_CHECK_MATRIX(k->x_);
    KALMAN_NULL_CHECK_MATRIX(k->P_);
    KALMAN_NULL_CHECK_MATRIX(k->H_);
    KALMAN_NULL_CHECK_MATRIX(k->R_);
    KALMAN_NULL_CHECK_MATRIX(k->z_);
    KALMAN_NULL_CHECK_MATRIX(k->H_TRANSPOSE);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_X_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FP_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_FPFT_);
    KALMAN_NULL_CHECK_MATRIX(k->TEMP_PHT_);

    // PHT = P * H^T
    KALMAN_RETURN_CHECK(transposeMatrix(k->H_, k->H_TRANSPOSE));
    KALMAN_RETURN_CHECK(multMatrix(k->P_, k->H_TRANSPOSE, k->TEMP_PHT_));

    // S = H * PHT + R
    KALMAN_RETURN_CHECK(multMatrix(k->H_, k->TEMP_PHT_, k->TEMP_FPFT_));
    KALMAN_RETURN_CHECK(addMatrix(k->TEMP_FPFT_, k->R_, k->TEMP_FP_));

    // S^-1
    matrixReturnCodes invRet =
        inverseMatrixWithJitter(k->TEMP_FP_, k->TEMP_FPFT_, (matrixType) 1e-6, 3, (matrixType) 100);
    if (invRet != MATRIX_SUCCESS)
    {
        LOG_ERROR("kalmanFilterUpdate failed to invert innovation matrix S.");
        return invRet;
    }

    // K = PHT * S^-1
    KALMAN_RETURN_CHECK(multMatrix(k->TEMP_PHT_, k->TEMP_FPFT_, k->TEMP_FP_));

    // residual = z - Hx (store back into z_)
    KALMAN_RETURN_CHECK(multMatrix(k->H_, k->x_, k->TEMP_X_));
    int i;
    for (i = 0; i < k->z_->row; ++i)
    {
        matrixType residual = ACCESS_MATRIX(*k->z_, i, 0) - ACCESS_MATRIX(*k->TEMP_X_, i, 0);
        SET_MATRIX(*k->z_, i, 0, residual);
    }

    // x = x + K * residual
    KALMAN_RETURN_CHECK(multMatrix(k->TEMP_FP_, k->z_, k->TEMP_X_));
    for (i = 0; i < k->x_->row; ++i)
    {
        matrixType updated = ACCESS_MATRIX(*k->x_, i, 0) + ACCESS_MATRIX(*k->TEMP_X_, i, 0);
        SET_MATRIX(*k->x_, i, 0, updated);
    }

    // P = (I - K H) P
    KALMAN_RETURN_CHECK(multMatrix(k->TEMP_FP_, k->H_, k->TEMP_PHT_));
    KALMAN_RETURN_CHECK(setIdentityMatrix(k->TEMP_FPFT_));
    // Store (I - K H) in TEMP_FP_ (K no longer needed after state update).
    KALMAN_RETURN_CHECK(subMatrix(k->TEMP_FPFT_, k->TEMP_PHT_, k->TEMP_FP_));
    KALMAN_RETURN_CHECK(multMatrix(k->TEMP_FP_, k->P_, k->TEMP_PHT_));
    KALMAN_RETURN_CHECK(copyMatrix(k->TEMP_PHT_, k->P_));

    return MATRIX_SUCCESS;
}
