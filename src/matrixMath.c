//#define TEST_MULT

#include "matrixMath.h"

// isnan check
#include <math.h>

void gaussianElimination(struct matrix* a, struct matrix* idenity, struct matrix* res);

#define MATRIX_INVERSE_STACK_MAX_N 16

// Inputs are three matrixes. a and b matrix will be multiplied and the result will be output to res
matrixReturnCodes multMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    int i, j, k;

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

    ZEROIZE_MATRIX(*res, i, j);

    // Do the multiplication
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            for (k = 0; k < b->col; ++k)
            {
                SET_MATRIX(*res, i, k, ACCESS_MATRIX(*res, i, k) + ACCESS_MATRIX(*a, i, j) * ACCESS_MATRIX(*b, j, k));
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes scaleMatrix(struct matrix* a, struct matrix* res, matrixType scaler)
{
    int i, j;

    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX_RES(res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    DIMENSION_CHECK_SCALER_MATRIX(a, res);

    NAN_CHECK_MATRIX(a);

    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            SET_MATRIX(*res, i, j, ACCESS_MATRIX(*a, i, j) * scaler);
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes addMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    int arow, acol;

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

    for (arow = 0; arow < a->row; ++arow)
    {
        for (acol = 0; acol < a->col; ++acol)
        {
            SET_MATRIX(*res, arow, acol, ACCESS_MATRIX(*a, arow, acol) + ACCESS_MATRIX(*b, arow, acol));
        }
    }

    return MATRIX_SUCCESS;
}

void printMatrix(struct matrix* a)
{
    int i, j;

    LOG_FUNCTION();

    if (a == NULL || !a->initilized)
    {
        printf("Trying to print a uninitilized matrix");
        return;
    }
    for (i = 0; i < a->row; ++i)
    {
        if (i > 0)
            printf("\n");
        for (j = 0; j < a->col; ++j)
        {
            // Cast to double so printing is safe for float/double/int matrixType overrides.
            printf("%f", (double) ACCESS_MATRIX(*a, i, j));
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
    int i, j;

    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NULL_CHECK_MATRIX(b);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(b);
    if (a->row != b->row || a->col != b->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
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
    int i, j;

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
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            SET_MATRIX(*b, j, i, ACCESS_MATRIX(*a, i, j));
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes subMatrix(struct matrix* a, struct matrix* b, struct matrix* res)
{
    int arow, acol;

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

    for (arow = 0; arow < a->row; ++arow)
    {
        for (acol = 0; acol < a->col; ++acol)
        {
            SET_MATRIX(*res, arow, acol, ACCESS_MATRIX(*a, arow, acol) - ACCESS_MATRIX(*b, arow, acol));
        }
    }

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrix(struct matrix* a, struct matrix* res)
{
    int i, j;

    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NULL_CHECK_MATRIX_RES(res);
    NO_ALIAS_CHECK_MATRIX2(a, res);

    NON_INIT_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(res);

    NAN_CHECK_MATRIX(a);

    // Zero out the result matrix
    for (i = 0; i < res->row; ++i)
    {
        for (j = 0; j < res->col; ++j)
        {
            SET_MATRIX(*res, i, j, 0);
        }
    }

    // Inverse for square matrix case
    if (a->row == a->col)
    {
        int n = a->row;

        /*
         * Fast path for small matrices: avoid heap alloc/free in the hot EKF update loop.
         * EKF innovation S is 5x5 in this project, so this path is always used there.
         */
        if (n <= MATRIX_INVERSE_STACK_MAX_N)
        {
            matrixType tempStorage[MATRIX_INVERSE_STACK_MAX_N * MATRIX_INVERSE_STACK_MAX_N];
            struct matrix tempMatrix;
            int r;
            int c;

            tempMatrix.mat = NULL;
            tempMatrix._mat = tempStorage;
            tempMatrix.initilized = 1;
            tempMatrix.row = n;
            tempMatrix.col = n;
            tempMatrix.jaggedAlloc = false;

            for (r = 0; r < n; ++r)
            {
                for (c = 0; c < n; ++c)
                {
                    SET_MATRIX(tempMatrix, r, c, 0);
                }
            }

            gaussianElimination(a, &tempMatrix, res);
        }
        else
        {
            // Fallback for larger matrices.
            struct matrix* tempMatrix = NULL;
            INIT_MATRIX(tempMatrix, a->row, a->col);
            gaussianElimination(a, tempMatrix, res);
            FREE_MATRIX(tempMatrix);
        }
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
    matrixReturnCodes invRet = MATRIX_ERROR;
    int attempt;
    int i;

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

    for (attempt = 0; attempt < maxAttempts; ++attempt)
    {
        invRet = inverseMatrix(a, res);
        if (invRet == MATRIX_SUCCESS)
        {
            return MATRIX_SUCCESS;
        }
        for (i = 0; i < a->row; ++i)
        {
            SET_MATRIX(*a, i, i, ACCESS_MATRIX(*a, i, i) + jitter);
        }
        jitter *= jitterScale;
    }
    return invRet;
}

matrixReturnCodes setIdentityMatrix(struct matrix* a)
{
    int i;
    int j;

    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);
    NON_INIT_CHECK_MATRIX(a);

    if (a->row != a->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }

    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            matrixType val = (i == j) ? (matrixType) 1 : (matrixType) 0;
            SET_MATRIX(*a, i, j, val);
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
                matrixType val = (i == j) ? (1 - ACCESS_MATRIX(*a, i, j)) : (-ACCESS_MATRIX(*a, i, j));
                SET_MATRIX(*res, i, j, val);
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
    int i, j;

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
    for (i = 0; i < a->row; ++i)
    {
        for (j = 0; j < a->col; ++j)
        {
            SET_MATRIX(*res, i, j, ACCESS_MATRIX(*a, i, j));
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes nanCheckMatrix(struct matrix* a)
{
    int i, j;

    LOG_FUNCTION();

    NULL_CHECK_MATRIX(a);

    NON_INIT_CHECK_MATRIX(a);

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
    matrixType temp;
    matrixType ratio;
    int i;
    int j;
    int k;

    LOG_FUNCTION();

    // Copy matrix a to res
    copyMatrix(a, res);

    // Initialize the inverse matrix as an identity matrix
    for (i = 0; i < res->col; ++i)
    {
        for (j = 0; j < res->col; ++j)
        {
            SET_MATRIX(*idenity, i, j, (i == j) ? (matrixType) 1 : (matrixType) 0);
        }
    }

    // Gaussian elimination
    for (i = 0; i < a->row; ++i)
    {
        temp = ACCESS_MATRIX(*res, i, i);
        for (j = 0; j < a->row; ++j)
        {
            SET_MATRIX(*res, i, j, ACCESS_MATRIX(*res, i, j) / temp);
            SET_MATRIX(*idenity, i, j, ACCESS_MATRIX(*idenity, i, j) / temp);
        }

        // Subtract the current row from all the other rows
        for (k = 0; k < a->row; ++k)
        {
            if (k != i)
            {
                ratio = ACCESS_MATRIX(*res, k, i);
                for (j = 0; j < a->row; ++j)
                {
                    SET_MATRIX(*res, k, j, ACCESS_MATRIX(*res, k, j) - ratio * ACCESS_MATRIX(*res, i, j));
                    SET_MATRIX(*idenity, k, j, ACCESS_MATRIX(*idenity, k, j) - ratio * ACCESS_MATRIX(*idenity, i, j));
                }
            }
        }
    }

    // Copy the inverse matrix to res
    copyMatrix(idenity, res);
}
