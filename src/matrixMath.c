//#define TEST_MULT

#include "matrixMath.h"

// isnan check
#include <math.h>

void gaussianElimination(struct matrix* a, struct matrix* idenity, struct matrix* res);

// Inputs are three matrixes. a and b matrix will be multiplied and the result will be output to res
matrixReturnCodes multMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX3(a, b, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);
    NON_INIT_CHECK_MATRIX(res);

    DIMENSION_CHECK_MULT_MATRIX(a, b, res);

    NAN_CHECK_MATRIX(a);
    NAN_CHECK_MATRIX(b);

    int i, j, k;
    ZEROIZE_MATRIX(*res, i, j);

    // Do the multiplication
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            for (k = 0; k < b->col; ++k)
            {
                if (res->jaggedAlloc)
                {
                    res->mat[i][k] += ACCESS_MATRIX(*a, i, j) * ACCESS_MATRIX(*b, j, k);
                }
                else
                {
                    ACCESS_STATIC_MATRIX(*res, i, k) += ACCESS_MATRIX(*a, i, j) * ACCESS_MATRIX(*b, j, k);
                }
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes scaleMatrix(struct matrix* a, struct matrix* res, matrixType scaler)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX_RES(res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    DIMENSION_CHECK_SCALER_MATRIX(a, res);

    NAN_CHECK_MATRIX(a);

    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            if (res->jaggedAlloc)
            {
                res->mat[i][j] = ACCESS_MATRIX(*a, i, j) * scaler;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, i, j) = ACCESS_MATRIX(*a, i, j) * scaler;
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes addMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX3(a, b, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);
    NON_INIT_CHECK_MATRIX(res);

    DIMENSION_CHECK_ADD_MATRIX(a, b, res);

    NAN_CHECK_MATRIX(a);
    NAN_CHECK_MATRIX(b);

    int arow, acol;
    for (arow = 0; arow < a->row; ++arow)
    {
        for (acol = 0; acol < a->col; ++acol)
        {
            if (res->jaggedAlloc)
            {
                res->mat[arow][acol] = ACCESS_MATRIX(*a, arow, acol) + ACCESS_MATRIX(*b, arow, acol);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, arow, acol) = ACCESS_MATRIX(*a, arow, acol) + ACCESS_MATRIX(*b, arow, acol);
            }
        }
    }

    return MATRIX_SUCCESS;
}

void printMatrix(struct matrix* a)
{
    LOG_FUNCTION();

    if (a == NULL || !a->initilized)
    {
        printf("Trying to print a uninitilized matrix");
        return;
    }
    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        if (i > 0)
            printf("\n");
        for (j = 0; j < a->col; ++j)
        {
#if defined(_INT)
            printf("%d", ACCESS_MATRIX(*a, i, j));
#elif defined(_DOUBLE)
            printf("%f", ACCESS_MATRIX(*a, i, j));
#else
            printf("           Type not defined!");
            printf("%d", ACCESS_MATRIX(*a, i, j));
#endif
            if (j < a->col - 1)
            {
                printf(",\t");
            }
        }
    }
    printf("\n");
}

matrixReturnCodes compareMatrieces(struct matrix* a, struct matrix* b)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);
    if (a->row != b->row || a->col != b->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            if (ACCESS_MATRIX(*a, i, j) != ACCESS_MATRIX(*b, i, j))
            {
                return MATRIX_COMPARE_FAILURE;
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes transposeMatrix(struct matrix* a, struct matrix* b)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);
    NO_ALIAS_CHECK_MATRIX2(a, b);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);

    NAN_CHECK_MATRIX(a);

    if (a->row != b->col || a->col != b->row)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            if (b->jaggedAlloc)
            {
                b->mat[j][i] = ACCESS_MATRIX(*a, i, j);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*b, j, i) = ACCESS_MATRIX(*a, i, j);
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes subMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX3(a, b, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);
    NON_INIT_CHECK_MATRIX(res);

    DIMENSION_CHECK_ADD_MATRIX(a, b, res);

    NAN_CHECK_MATRIX(a);
    NAN_CHECK_MATRIX(b);

    int arow, acol;
    for (arow = 0; arow < a->row; ++arow)
    {
        for (acol = 0; acol < a->col; ++acol)
        {
            if (res->jaggedAlloc)
            {
                res->mat[arow][acol] = ACCESS_MATRIX(*a, arow, acol) - ACCESS_MATRIX(*b, arow, acol);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, arow, acol) = ACCESS_MATRIX(*a, arow, acol) - ACCESS_MATRIX(*b, arow, acol);
            }
        }
    }

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrix(struct matrix* a, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX2(a, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    NAN_CHECK_MATRIX(a);

    // Zero out the result matrix
    int i, j;
    for (i = 0; i < res->row; ++i)
    {
        for (j = 0; j < res->col; ++j)
        {
            if (res->jaggedAlloc)
            {
                res->mat[i][j] = 0;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, i, j) = 0;
            }
        }
    }

    // Inverse for square matrax case
    if (a->row == a->col)
    {
        // We need to make a temp matrix to store the identity matrix
        // TODO: Make this a static matrix that is in the EKF struct and passed in
        struct matrix* tempMatrix = NULL;
        INIT_MATRIX(tempMatrix, a->row, a->col);
        // Gaussian elimination
        gaussianElimination(a, tempMatrix, res);
        FREE_MATRIX(tempMatrix);
    }
    else
    {
        // TODO: Implement inverse for non square matrix
        return MATRIX_DIMENSION_MISMATCH;
    }

    NAN_CHECK_MATRIX(res);

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrixWithJitter(struct matrix* a, struct matrix* res, matrixType jitter, int maxAttempts,
                                          matrixType jitterScale)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX2(a, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    NAN_CHECK_MATRIX(a);

    if (a->row != a->col || res->row != res->col || a->row != res->row)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    if (maxAttempts <= 0 || jitter <= 0)
    {
        return MATRIX_ERROR;
    }

    matrixReturnCodes invRet = MATRIX_ERROR;
    for (int attempt = 0; attempt < maxAttempts; ++attempt)
    {
        invRet = inverseMatrix(a, res);
        if (invRet == MATRIX_SUCCESS)
        {
            return MATRIX_SUCCESS;
        }
        for (int i = 0; i < a->row; ++i)
        {
            if (a->jaggedAlloc)
            {
                a->mat[i][i] += jitter;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*a, i, i) += jitter;
            }
        }
        jitter *= jitterScale;
    }
    return invRet;
}

matrixReturnCodes setIdentityMatrix(struct matrix* a)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(a);

    if (a->row != a->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }

    for (int i = 0; i < a->row; ++i)
    {
        for (int j = 0; j < a->col; ++j)
        {
            matrixType val = (i == j) ? 1 : 0;
            if (a->jaggedAlloc)
            {
                a->mat[i][j] = val;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*a, i, j) = val;
            }
        }
    }

    return MATRIX_SUCCESS;
}

matrixReturnCodes identityMatrixMinusA(struct matrix* a, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX2(a, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    NAN_CHECK_MATRIX(a);

    // Inverse for square matrax case
    if (a->row == a->col)
    {
        // Since we dont want to allocate a whole another matrix we will reimplemnt the subtract function
        int i, j;
        for (i = 0; i < a->row; ++i)
        {
            for (j = 0; j < a->col; ++j)
            {
                if (res->jaggedAlloc)
                {
                    if (i == j)
                    {
                        res->mat[i][j] = 1 - ACCESS_MATRIX(*a, i, j);
                    }
                    else
                    {
                        res->mat[i][j] = -ACCESS_MATRIX(*a, i, j);
                    }
                }
                else
                {
                    if (i == j)
                    {
                        ACCESS_STATIC_MATRIX(*res, i, j) = 1 - ACCESS_MATRIX(*a, i, j);
                    }
                    else
                    {
                        ACCESS_STATIC_MATRIX(*res, i, j) = -ACCESS_MATRIX(*a, i, j);
                    }
                }
            }
        }
    }
    else
    {
        // TODO not sure if this is needed
        return MATRIX_DIMENSION_MISMATCH;
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes copyMatrix(struct matrix* a, struct matrix* res)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX2(a, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    if (a->row != res->row || a->col != res->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            if (res->jaggedAlloc)
            {
                res->mat[i][j] = ACCESS_MATRIX(*a, i, j);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, i, j) = ACCESS_MATRIX(*a, i, j);
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes nanCheckMatrix(struct matrix* a)
{
    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NON_INIT_CHECK_MATRIX(a);

    int i, j;
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            if (isnan(ACCESS_MATRIX(*a, i, j)))
            {
                return MATRIX_NAN_FAILURE;
            }
            if (isinf(ACCESS_MATRIX(*a, i, j)))
            {
                return MATRIX_INF_FAILURE;
            }
        }
    }
    return MATRIX_SUCCESS;
}

void gaussianElimination(struct matrix* a, struct matrix* idenity, struct matrix* res)
{
    LOG_FUNCTION();

    matrixType temp;
    matrixType ratio;

    // Copy matrix a to res
    copyMatrix(a, res);

    // Initialize the inverse matrix as an identity matrix
    for (int i = 0; i < res->col; ++i)
    {
        for (int j = 0; j < res->col; ++j)
        {
            if (i == j)
            {
                if (idenity->jaggedAlloc)
                {
                    idenity->mat[i][j] = 1;
                }
                else
                {
                    ACCESS_STATIC_MATRIX(*idenity, i, j) = 1;
                }
            }
            else
            {
                if (idenity->jaggedAlloc)
                {
                    idenity->mat[i][j] = 0;
                }
                else
                {
                    ACCESS_STATIC_MATRIX(*idenity, i, j) = 0;
                }
            }
        }
    }

    // Gaussian elimination
    for (int i = 0; i < a->row; ++i)
    {
        temp = ACCESS_MATRIX(*res, i, i);
        for (int j = 0; j < a->row; ++j)
        {
            if (res->jaggedAlloc)
            {
                res->mat[i][j] /= temp;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*res, i, j) /= temp;
            }
            if (idenity->jaggedAlloc)
            {
                idenity->mat[i][j] /= temp;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*idenity, i, j) /= temp;
            }
        }

        // Subtract the current row from all the other rows
        for (int k = 0; k < a->row; ++k)
        {
            if (k != i)
            {
                ratio = ACCESS_MATRIX(*res, k, i);
                for (int j = 0; j < a->row; ++j)
                {
                    if (res->jaggedAlloc)
                    {
                        res->mat[k][j] -= ratio * ACCESS_MATRIX(*res, i, j);
                    }
                    else
                    {
                        ACCESS_STATIC_MATRIX(*res, k, j) -= ratio * ACCESS_MATRIX(*res, i, j);
                    }
                    if (idenity->jaggedAlloc)
                    {
                        idenity->mat[k][j] -= ratio * ACCESS_MATRIX(*idenity, i, j);
                    }
                    else
                    {
                        ACCESS_STATIC_MATRIX(*idenity, k, j) -= ratio * ACCESS_MATRIX(*idenity, i, j);
                    }
                }
            }
        }
    }

    // Copy the inverse matrix to res
    copyMatrix(idenity, res);
}
