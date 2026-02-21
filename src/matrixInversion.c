#include "matrixMath.h"

#include <math.h>
#include <stdlib.h>

static matrixType matrixInversionEpsilon(void)
{
    return (matrixType) 1e-12;
}

static matrixType matrixAbs(matrixType value)
{
    return (matrixType) fabs((double) value);
}

static matrixReturnCodes validateInverseInputs(struct matrix* a, struct matrix* res)
{
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

    return MATRIX_SUCCESS;
}

static void swapRows(struct matrix* m, int rowA, int rowB)
{
    int j;

    if (rowA == rowB)
    {
        return;
    }

    for (j = 0; j < m->col; ++j)
    {
        matrixType tmp = ACCESS_MATRIX(*m, rowA, j);
        SET_MATRIX(*m, rowA, j, ACCESS_MATRIX(*m, rowB, j));
        SET_MATRIX(*m, rowB, j, tmp);
    }
}

static bool isSymmetricWithinTolerance(struct matrix* a, matrixType tol)
{
    int i;
    int j;

    for (i = 0; i < a->row; ++i)
    {
        for (j = i + 1; j < a->col; ++j)
        {
            matrixType diff = matrixAbs(ACCESS_MATRIX(*a, i, j) - ACCESS_MATRIX(*a, j, i));
            if (diff > tol)
            {
                return false;
            }
        }
    }

    return true;
}

static matrixReturnCodes inverseMatrixGaussJordanImpl(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;
    matrixType eps = matrixInversionEpsilon();
    int n = a->row;
    int i;
    int r;
    int c;
    struct matrix* work = NULL;

    INIT_MATRIX(work, n, n);

    ret = copyMatrix(a, work);
    if (ret != MATRIX_SUCCESS)
    {
        FREE_MATRIX(work);
        return ret;
    }

    ret = setIdentityMatrix(res);
    if (ret != MATRIX_SUCCESS)
    {
        FREE_MATRIX(work);
        return ret;
    }

    for (i = 0; i < n; ++i)
    {
        int pivotRow = i;
        matrixType pivotAbs = matrixAbs(ACCESS_MATRIX(*work, i, i));

        for (r = i + 1; r < n; ++r)
        {
            matrixType candidateAbs = matrixAbs(ACCESS_MATRIX(*work, r, i));
            if (candidateAbs > pivotAbs)
            {
                pivotAbs = candidateAbs;
                pivotRow = r;
            }
        }

        if (pivotAbs <= eps)
        {
            FREE_MATRIX(work);
            return MATRIX_ERROR;
        }

        swapRows(work, i, pivotRow);
        swapRows(res, i, pivotRow);

        {
            matrixType pivot = ACCESS_MATRIX(*work, i, i);
            for (c = 0; c < n; ++c)
            {
                SET_MATRIX(*work, i, c, ACCESS_MATRIX(*work, i, c) / pivot);
                SET_MATRIX(*res, i, c, ACCESS_MATRIX(*res, i, c) / pivot);
            }
        }

        for (r = 0; r < n; ++r)
        {
            if (r != i)
            {
                matrixType factor = ACCESS_MATRIX(*work, r, i);
                if (matrixAbs(factor) > eps)
                {
                    for (c = 0; c < n; ++c)
                    {
                        SET_MATRIX(*work, r, c,
                                   ACCESS_MATRIX(*work, r, c) - factor * ACCESS_MATRIX(*work, i, c));
                        SET_MATRIX(*res, r, c, ACCESS_MATRIX(*res, r, c) - factor * ACCESS_MATRIX(*res, i, c));
                    }
                }
                SET_MATRIX(*work, r, i, 0);
            }
        }
    }

    FREE_MATRIX(work);
    return MATRIX_SUCCESS;
}

static matrixReturnCodes inverseMatrixLUImpl(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;
    matrixType eps = matrixInversionEpsilon();
    int n = a->row;
    int i;
    int j;
    int k;
    int col;
    int* perm = NULL;
    matrixType* rhs = NULL;
    matrixType* y = NULL;
    matrixType* x = NULL;
    struct matrix* lu = NULL;

    INIT_MATRIX(lu, n, n);

    ret = copyMatrix(a, lu);
    if (ret != MATRIX_SUCCESS)
    {
        FREE_MATRIX(lu);
        return ret;
    }

    perm = malloc((size_t) n * sizeof(int));
    rhs = malloc((size_t) n * sizeof(matrixType));
    y = malloc((size_t) n * sizeof(matrixType));
    x = malloc((size_t) n * sizeof(matrixType));

    if (perm == NULL || rhs == NULL || y == NULL || x == NULL)
    {
        FREE_MATRIX(lu);
        free(perm);
        free(rhs);
        free(y);
        free(x);
        return MATRIX_ERROR;
    }

    for (i = 0; i < n; ++i)
    {
        perm[i] = i;
    }

    for (k = 0; k < n; ++k)
    {
        int pivotRow = k;
        matrixType pivotAbs = matrixAbs(ACCESS_MATRIX(*lu, k, k));

        for (i = k + 1; i < n; ++i)
        {
            matrixType candidateAbs = matrixAbs(ACCESS_MATRIX(*lu, i, k));
            if (candidateAbs > pivotAbs)
            {
                pivotAbs = candidateAbs;
                pivotRow = i;
            }
        }

        if (pivotAbs <= eps)
        {
            FREE_MATRIX(lu);
            free(perm);
            free(rhs);
            free(y);
            free(x);
            return MATRIX_ERROR;
        }

        if (pivotRow != k)
        {
            int tmpPerm;
            swapRows(lu, k, pivotRow);
            tmpPerm = perm[k];
            perm[k] = perm[pivotRow];
            perm[pivotRow] = tmpPerm;
        }

        for (i = k + 1; i < n; ++i)
        {
            matrixType pivot = ACCESS_MATRIX(*lu, k, k);
            matrixType multiplier = ACCESS_MATRIX(*lu, i, k) / pivot;
            SET_MATRIX(*lu, i, k, multiplier);

            for (j = k + 1; j < n; ++j)
            {
                SET_MATRIX(*lu, i, j,
                           ACCESS_MATRIX(*lu, i, j) - multiplier * ACCESS_MATRIX(*lu, k, j));
            }
        }
    }

    for (col = 0; col < n; ++col)
    {
        for (i = 0; i < n; ++i)
        {
            rhs[i] = (perm[i] == col) ? (matrixType) 1 : (matrixType) 0;
        }

        for (i = 0; i < n; ++i)
        {
            matrixType sum = rhs[i];
            for (j = 0; j < i; ++j)
            {
                sum -= ACCESS_MATRIX(*lu, i, j) * y[j];
            }
            y[i] = sum;
        }

        for (i = n - 1; i >= 0; --i)
        {
            matrixType sum = y[i];
            matrixType pivot = ACCESS_MATRIX(*lu, i, i);

            for (j = i + 1; j < n; ++j)
            {
                sum -= ACCESS_MATRIX(*lu, i, j) * x[j];
            }

            if (matrixAbs(pivot) <= eps)
            {
                FREE_MATRIX(lu);
                free(perm);
                free(rhs);
                free(y);
                free(x);
                return MATRIX_ERROR;
            }

            x[i] = sum / pivot;
        }

        for (i = 0; i < n; ++i)
        {
            SET_MATRIX(*res, i, col, x[i]);
        }
    }

    FREE_MATRIX(lu);
    free(perm);
    free(rhs);
    free(y);
    free(x);

    return MATRIX_SUCCESS;
}

static matrixReturnCodes inverseMatrixCholeskyImpl(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;
    matrixType eps = matrixInversionEpsilon();
    matrixType symmetryTol = (matrixType) (10 * eps);
    int n = a->row;
    int i;
    int j;
    int k;
    int col;
    matrixType* y = NULL;
    matrixType* x = NULL;
    struct matrix* L = NULL;

    if (!isSymmetricWithinTolerance(a, symmetryTol))
    {
        return MATRIX_ERROR;
    }

    INIT_MATRIX(L, n, n);

    y = malloc((size_t) n * sizeof(matrixType));
    x = malloc((size_t) n * sizeof(matrixType));
    if (y == NULL || x == NULL)
    {
        FREE_MATRIX(L);
        free(y);
        free(x);
        return MATRIX_ERROR;
    }

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j <= i; ++j)
        {
            matrixType sum = ACCESS_MATRIX(*a, i, j);
            for (k = 0; k < j; ++k)
            {
                sum -= ACCESS_MATRIX(*L, i, k) * ACCESS_MATRIX(*L, j, k);
            }

            if (i == j)
            {
                if (sum <= eps)
                {
                    FREE_MATRIX(L);
                    free(y);
                    free(x);
                    return MATRIX_ERROR;
                }
                SET_MATRIX(*L, i, j, (matrixType) sqrt((double) sum));
            }
            else
            {
                matrixType diag = ACCESS_MATRIX(*L, j, j);
                if (matrixAbs(diag) <= eps)
                {
                    FREE_MATRIX(L);
                    free(y);
                    free(x);
                    return MATRIX_ERROR;
                }
                SET_MATRIX(*L, i, j, sum / diag);
            }
        }
    }

    for (col = 0; col < n; ++col)
    {
        for (i = 0; i < n; ++i)
        {
            matrixType rhs = (i == col) ? (matrixType) 1 : (matrixType) 0;
            matrixType sum = rhs;
            matrixType diag = ACCESS_MATRIX(*L, i, i);

            for (k = 0; k < i; ++k)
            {
                sum -= ACCESS_MATRIX(*L, i, k) * y[k];
            }

            if (matrixAbs(diag) <= eps)
            {
                FREE_MATRIX(L);
                free(y);
                free(x);
                return MATRIX_ERROR;
            }

            y[i] = sum / diag;
        }

        for (i = n - 1; i >= 0; --i)
        {
            matrixType sum = y[i];
            matrixType diag = ACCESS_MATRIX(*L, i, i);

            for (k = i + 1; k < n; ++k)
            {
                sum -= ACCESS_MATRIX(*L, k, i) * x[k];
            }

            if (matrixAbs(diag) <= eps)
            {
                FREE_MATRIX(L);
                free(y);
                free(x);
                return MATRIX_ERROR;
            }

            x[i] = sum / diag;
        }

        for (i = 0; i < n; ++i)
        {
            SET_MATRIX(*res, i, col, x[i]);
        }
    }

    for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            matrixType avg = (matrixType) (0.5 * (ACCESS_MATRIX(*res, i, j) + ACCESS_MATRIX(*res, j, i)));
            SET_MATRIX(*res, i, j, avg);
            SET_MATRIX(*res, j, i, avg);
        }
    }

    FREE_MATRIX(L);
    free(y);
    free(x);

    ret = nanCheckMatrix(res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrixGaussJordan(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;

    LOG_FUNCTION();

    ret = validateInverseInputs(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    ret = inverseMatrixGaussJordanImpl(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    NAN_CHECK_MATRIX(res);

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrixLU(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;

    LOG_FUNCTION();

    ret = validateInverseInputs(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    ret = inverseMatrixLUImpl(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    NAN_CHECK_MATRIX(res);

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrixCholesky(struct matrix* a, struct matrix* res)
{
    matrixReturnCodes ret;

    LOG_FUNCTION();

    ret = validateInverseInputs(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    ret = inverseMatrixCholeskyImpl(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    NAN_CHECK_MATRIX(res);

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrixByMethod(struct matrix* a, struct matrix* res, matrixInversionMethod method)
{
    matrixReturnCodes ret;

    LOG_FUNCTION();

    ret = validateInverseInputs(a, res);
    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    switch (method)
    {
        case MATRIX_INVERSE_GAUSS_JORDAN:
            ret = inverseMatrixGaussJordanImpl(a, res);
            break;

        case MATRIX_INVERSE_LU:
            ret = inverseMatrixLUImpl(a, res);
            break;

        case MATRIX_INVERSE_CHOLESKY:
            ret = inverseMatrixCholeskyImpl(a, res);
            break;

        case MATRIX_INVERSE_AUTO:
            ret = inverseMatrixCholeskyImpl(a, res);
            if (ret == MATRIX_SUCCESS)
            {
                break;
            }
            ret = inverseMatrixLUImpl(a, res);
            if (ret == MATRIX_SUCCESS)
            {
                break;
            }
            ret = inverseMatrixGaussJordanImpl(a, res);
            break;

        default:
            return MATRIX_ERROR;
    }

    if (ret != MATRIX_SUCCESS)
    {
        return ret;
    }

    NAN_CHECK_MATRIX(res);

    return MATRIX_SUCCESS;
}
