/**
 * @file ekfMatrixMathTest.c
 * 
 * I know this may leak memory if malloc fails and all other things... Its a test. I might fix it I might not.
*/

#include "EKF.h"
#include <math.h>
#include <string.h>

// Tests return matrixReturnCodes; override EKF macro to avoid enum-conversion warnings.
#ifdef MATRIX_MATH_RETURN_CHECK
#undef MATRIX_MATH_RETURN_CHECK
#endif
#define MATRIX_MATH_RETURN_CHECK(ret) \
    do { \
        matrixReturnCodes _ret = (ret); \
        if (_ret != MATRIX_SUCCESS) { \
            return _ret; \
        } \
    } while (0)

matrixReturnCodes copyTest(bool increasedLogging);
matrixReturnCodes inverseTest(bool increasedLogging);
matrixReturnCodes inverseTest3x3(bool increasedLogging);
matrixReturnCodes inverseTest4x4(bool increasedLogging);
matrixReturnCodes inverseTest5x5(bool increasedLogging);
matrixReturnCodes inverseTest6x6(bool increasedLogging);
matrixReturnCodes additionTest(bool increasedLogging);
matrixReturnCodes subtractionTest(bool increasedLogging);
matrixReturnCodes multiplicationTest(bool increasedLogging);
matrixReturnCodes tranposeTest(bool increasedLogging);
matrixReturnCodes scaleTest(bool increasedLogging);
matrixReturnCodes idenityMatrixMinusATest(bool increasedLogging);
matrixReturnCodes dimensionMismatchTests(bool increasedLogging);
matrixReturnCodes nanAndInfTests(bool increasedLogging);
matrixReturnCodes staticMatrixAddTest(bool increasedLogging);
matrixReturnCodes staticMatrixMultTest(bool increasedLogging);
matrixReturnCodes staticMatrixTransposeTest(bool increasedLogging);
matrixReturnCodes staticMatrixScaleTest(bool increasedLogging);
matrixReturnCodes staticMatrixIdentityMinusATest(bool increasedLogging);
matrixReturnCodes compareDimensionMismatchTest(bool increasedLogging);
matrixReturnCodes copy2DArrayTest(bool increasedLogging);
matrixReturnCodes multMatrixZeroingTest(bool increasedLogging);
matrixReturnCodes inverseSingularTest(bool increasedLogging);
matrixReturnCodes inverseJitterTest(bool increasedLogging);
matrixReturnCodes nanPropagationTests(bool increasedLogging);
matrixReturnCodes identityPropertiesTest(bool increasedLogging);
matrixReturnCodes setIdentityMatrixTest(bool increasedLogging);
matrixReturnCodes transposeTwiceTest(bool increasedLogging);
matrixReturnCodes scaleIdentityZeroTest(bool increasedLogging);
matrixReturnCodes inverseSanityTest(bool increasedLogging);
matrixReturnCodes inverseDoesNotMutateTest(bool increasedLogging);
matrixReturnCodes jaggedStaticParityTest(bool increasedLogging);
matrixReturnCodes nonsquareMultiplyTest(bool increasedLogging);
matrixReturnCodes nanIdentityMinusATest(bool increasedLogging);
matrixReturnCodes aliasingGuardTests(bool increasedLogging);

typedef struct TestCase_ {
    const char *name;
    matrixReturnCodes (*fn)(bool increasedLogging);
} TestCase;

static bool isMatch(const char *target, const char *arg)
{
    return (arg != NULL) && (target != NULL) && (strcmp(target, arg) == 0);
}

int main(int argc, char **argv)
{
    LOG_INFO("Starting Matrix Math Tests");
    bool increasedLogging = true;
    const char *filter = NULL;
    if (argc > 1)
    {
        filter = argv[1];
    }

    TestCase tests[] = {
        {"copy", copyTest},
        {"inverse2x2", inverseTest},
        {"inverse3x3", inverseTest3x3},
        {"inverse4x4", inverseTest4x4},
        {"inverse5x5", inverseTest5x5},
        {"inverse6x6", inverseTest6x6},
        {"add", additionTest},
        {"sub", subtractionTest},
        {"mult", multiplicationTest},
        {"transpose", tranposeTest},
        {"scale", scaleTest},
        {"identity_minus_a", idenityMatrixMinusATest},
        {"dimension_mismatch", dimensionMismatchTests},
        {"nan_inf_check", nanAndInfTests},
        {"static_add", staticMatrixAddTest},
        {"static_mult", staticMatrixMultTest},
        {"static_transpose", staticMatrixTransposeTest},
        {"static_scale", staticMatrixScaleTest},
        {"static_identity_minus_a", staticMatrixIdentityMinusATest},
        {"compare_mismatch", compareDimensionMismatchTest},
        {"copy_2d_array", copy2DArrayTest},
        {"mult_zeroing", multMatrixZeroingTest},
        {"inverse_singular", inverseSingularTest},
        {"inverse_jitter", inverseJitterTest},
        {"nan_propagation", nanPropagationTests},
        {"identity_properties", identityPropertiesTest},
        {"set_identity", setIdentityMatrixTest},
        {"transpose_twice", transposeTwiceTest},
        {"scale_identity_zero", scaleIdentityZeroTest},
        {"inverse_sanity", inverseSanityTest},
        {"inverse_no_mutate", inverseDoesNotMutateTest},
        {"jagged_static_parity", jaggedStaticParityTest},
        {"nonsquare_mult", nonsquareMultiplyTest},
        {"nan_identity_minus_a", nanIdentityMinusATest},
        {"aliasing_guards", aliasingGuardTests},
    };

    int testCount = (int)(sizeof(tests) / sizeof(tests[0]));
    for (int i = 0; i < testCount; ++i)
    {
        if (filter == NULL || isMatch(tests[i].name, filter))
        {
            MATRIX_MATH_RETURN_CHECK(tests[i].fn(increasedLogging));
        }
    }

    LOG_INFO("Completed All Tests");
    return 0;
}

matrixReturnCodes tranposeTest(bool increasedLogging)
{
    struct matrix *tTest = NULL;
    struct matrix *correctResult = NULL;
    struct matrix *result = NULL;

    INIT_MATRIX(tTest, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);
    INIT_MATRIX(result, 2, 2);

    LOG_INFO("Starting Transpose Test");
    tTest->mat[0][0] = 1;
    tTest->mat[0][1] = 2;
    tTest->mat[1][0] = 3;
    tTest->mat[1][1] = 4;
    correctResult->mat[0][0] = 1;
    correctResult->mat[0][1] = 3;
    correctResult->mat[1][0] = 2;
    correctResult->mat[1][1] = 4;

    MATRIX_MATH_RETURN_CHECK(transposeMatrix(tTest, result));
    if (increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(tTest);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    FREE_MATRIX(tTest);
    FREE_MATRIX(correctResult);
    FREE_MATRIX(result);

    LOG_INFO("Completed Transpose Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes scaleTest(bool increasedLogging)
{
    struct matrix *sTest = NULL;
    struct matrix *correctResult = NULL;
    struct matrix *result = NULL;

    INIT_MATRIX(sTest, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);
    INIT_MATRIX(result, 2, 2);

    LOG_INFO("Starting Scale Test");
    sTest->mat[0][0] = 1;
    sTest->mat[0][1] = 2;
    sTest->mat[1][0] = 3;
    sTest->mat[1][1] = 4;
    correctResult->mat[0][0] = 2;
    correctResult->mat[0][1] = 4;
    correctResult->mat[1][0] = 6;
    correctResult->mat[1][1] = 8;

    MATRIX_MATH_RETURN_CHECK(scaleMatrix(sTest, result, 2));
    if (increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(sTest);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    FREE_MATRIX(sTest);
    FREE_MATRIX(correctResult);
    FREE_MATRIX(result);

    LOG_INFO("Completed Scale Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes idenityMatrixMinusATest(bool increasedLogging)
{
    struct matrix *iTest = NULL;
    struct matrix *correctResult = NULL;
    struct matrix *result = NULL;

    INIT_MATRIX(iTest, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);
    INIT_MATRIX(result, 2, 2);

    LOG_INFO("Starting I - A Test");
    iTest->mat[0][0] = 1;
    iTest->mat[0][1] = 2;
    iTest->mat[1][0] = 3;
    iTest->mat[1][1] = 4;
    correctResult->mat[0][0] = 0;
    correctResult->mat[0][1] = -2;
    correctResult->mat[1][0] = -3;
    correctResult->mat[1][1] = -3;

    MATRIX_MATH_RETURN_CHECK(identityMatrixMinusA(iTest, result));
    if (increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(iTest);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    FREE_MATRIX(iTest);
    FREE_MATRIX(correctResult);
    FREE_MATRIX(result);

    LOG_INFO("Completed I - A Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes multiplicationTest(bool increasedLogging)
{
    struct matrix *multTest = NULL;
    struct matrix *multTest2 = NULL;

    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(multTest, 2, 2);
    INIT_MATRIX(multTest2, 2, 2);

    INIT_MATRIX(result, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);

    LOG_INFO("Starting Mult Test");
    multTest->mat[0][0] = 1;
    multTest->mat[0][1] = 2;
    multTest->mat[1][0] = 3;
    multTest->mat[1][1] = 4;
    multTest2->mat[0][0] = 5;
    multTest2->mat[0][1] = 6;
    multTest2->mat[1][0] = 7;
    multTest2->mat[1][1] = 8;
    correctResult->mat[0][0] = 19;
    correctResult->mat[0][1] = 22;
    correctResult->mat[1][0] = 43;
    correctResult->mat[1][1] = 50;

    MATRIX_MATH_RETURN_CHECK(multMatrix(multTest, multTest2, result));

    if(increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(multTest);
        LOG_INFO("Matrix 2");
        printMatrix(multTest2);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    FREE_MATRIX(correctResult);
    FREE_MATRIX(multTest);
    FREE_MATRIX(multTest2);
    FREE_MATRIX(result);

    return MATRIX_SUCCESS;
}

matrixReturnCodes subtractionTest(bool increasedLogging)
{
    struct matrix *subTest = NULL;
    struct matrix *subTest2 = NULL;
    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(subTest, 2, 2);
    INIT_MATRIX(subTest2, 2, 2);
    INIT_MATRIX(result, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);

    LOG_INFO("Starting Sub Test");
    subTest->mat[0][0] = 1;
    subTest->mat[0][1] = 2;
    subTest->mat[1][0] = 3;
    subTest->mat[1][1] = 4;
    subTest2->mat[0][0] = 5;
    subTest2->mat[0][1] = 6;
    subTest2->mat[1][0] = 7;
    subTest2->mat[1][1] = 8;
    correctResult->mat[0][0] = -4;
    correctResult->mat[0][1] = -4;
    correctResult->mat[1][0] = -4;
    correctResult->mat[1][1] = -4;

    MATRIX_MATH_RETURN_CHECK(subMatrix(subTest, subTest2, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(subTest);
        LOG_INFO("Matrix 2");
        printMatrix(subTest2);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(subTest);
    FREE_MATRIX(subTest2);
    FREE_MATRIX(result);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Sub Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes additionTest(bool increasedLogging)
{
    struct matrix *addTest = NULL;
    struct matrix *addTest2 = NULL;
    struct matrix *correctResult = NULL;
    struct matrix *result = NULL;

    INIT_MATRIX(addTest, 2, 2);
    INIT_MATRIX(addTest2, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);
    INIT_MATRIX(result, 2, 2);

    LOG_INFO("Starting Add Test");
    addTest->mat[0][0] = 1;
    addTest->mat[0][1] = 2;
    addTest->mat[1][0] = 3;
    addTest->mat[1][1] = 4;

    addTest2->mat[0][0] = 5;
    addTest2->mat[0][1] = 6;
    addTest2->mat[1][0] = 7;
    addTest2->mat[1][1] = 8;

    correctResult->mat[0][0] = 6;
    correctResult->mat[0][1] = 8;
    correctResult->mat[1][0] = 10;
    correctResult->mat[1][1] = 12;

    MATRIX_MATH_RETURN_CHECK(addMatrix(addTest, addTest2, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Matrix 1");
        printMatrix(addTest);
        LOG_INFO("Matrix 2");
        printMatrix(addTest2);
        LOG_INFO("Result");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(addTest);
    FREE_MATRIX(addTest2);
    FREE_MATRIX(correctResult);
    FREE_MATRIX(result);

    LOG_INFO("Completed Add Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseTest(bool increasedLogging)
{
    struct matrix *invTest = NULL;
    struct matrix *invTest2 = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(invTest, 2, 2);
    INIT_MATRIX(invTest2, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);

    LOG_INFO("Starting Inverse Test");
    invTest->mat[0][0] = 1;
    invTest->mat[0][1] = 2;
    invTest->mat[1][0] = 3;
    invTest->mat[1][1] = 4;

    correctResult->mat[0][0] = -2;
    correctResult->mat[0][1] = 1;
    correctResult->mat[1][0] = 1.5;
    correctResult->mat[1][1] = -0.5;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(invTest, invTest2));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(invTest2, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(invTest);
        LOG_INFO("Inverse");
        printMatrix(invTest2);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(invTest);
    FREE_MATRIX(invTest2);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Inverse Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseTest3x3(bool increasedLogging)
{
    struct matrix *invTest = NULL;
    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(invTest, 3, 3);
    INIT_MATRIX(result, 3, 3);
    INIT_MATRIX(correctResult, 3, 3);

    LOG_INFO("Starting Inverse 3x3 Test");
    invTest->mat[0][0] = 1;
    invTest->mat[0][1] = 2;
    invTest->mat[0][2] = 3;
    invTest->mat[1][0] = 0;
    invTest->mat[1][1] = 1;
    invTest->mat[1][2] = 4;
    invTest->mat[2][0] = 5;
    invTest->mat[2][1] = 6;
    invTest->mat[2][2] = 0;

    correctResult->mat[0][0] = -24;
    correctResult->mat[0][1] = 18;
    correctResult->mat[0][2] = 5;
    correctResult->mat[1][0] = 20;
    correctResult->mat[1][1] = -15;
    correctResult->mat[1][2] = -4;
    correctResult->mat[2][0] = -5;
    correctResult->mat[2][1] = 4;
    correctResult->mat[2][2] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(invTest, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(invTest);
        LOG_INFO("Inverse");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(invTest);
    FREE_MATRIX(result);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Inverse 3x3 Test");

    return MATRIX_SUCCESS;

}

matrixReturnCodes inverseTest4x4(bool increasedLogging)
{
    struct matrix *invTest = NULL;
    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(invTest, 4, 4);
    INIT_MATRIX(result, 4, 4);
    INIT_MATRIX(correctResult, 4, 4);

    LOG_INFO("Starting Inverse 4x4 Test");
    invTest->mat[0][0] = 1;
    invTest->mat[0][1] = 2;
    invTest->mat[0][2] = 3;
    invTest->mat[0][3] = 4;
    invTest->mat[1][0] = 0;
    invTest->mat[1][1] = 1;
    invTest->mat[1][2] = 2;
    invTest->mat[1][3] = 3;
    invTest->mat[2][0] = 0;
    invTest->mat[2][1] = 0;
    invTest->mat[2][2] = 1;
    invTest->mat[2][3] = 2;
    invTest->mat[3][0] = 0;
    invTest->mat[3][1] = 0;
    invTest->mat[3][2] = 0;
    invTest->mat[3][3] = 1;

    correctResult->mat[0][0] = 1;
    correctResult->mat[0][1] = -2;
    correctResult->mat[0][2] = 1;
    correctResult->mat[0][3] = 0;
    correctResult->mat[1][0] = 0;
    correctResult->mat[1][1] = 1;
    correctResult->mat[1][2] = -2;
    correctResult->mat[1][3] = 1;
    correctResult->mat[2][0] = 0;
    correctResult->mat[2][1] = 0;
    correctResult->mat[2][2] = 1;
    correctResult->mat[2][3] = -2;
    correctResult->mat[3][0] = 0;
    correctResult->mat[3][1] = 0;
    correctResult->mat[3][2] = 0;
    correctResult->mat[3][3] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(invTest, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(invTest);
        LOG_INFO("Inverse");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(invTest);
    FREE_MATRIX(result);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Inverse 4x4 Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseTest5x5(bool increasedLogging)
{
    struct matrix *invTest = NULL;
    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(invTest, 5, 5);
    INIT_MATRIX(result, 5, 5);
    INIT_MATRIX(correctResult, 5, 5);

    LOG_INFO("Starting Inverse 5x5 Test");
    invTest->mat[0][0] = 1;
    invTest->mat[0][1] = 2;
    invTest->mat[0][2] = 3;
    invTest->mat[0][3] = 4;
    invTest->mat[0][4] = 5;
    invTest->mat[1][0] = 0;
    invTest->mat[1][1] = 1;
    invTest->mat[1][2] = 2;
    invTest->mat[1][3] = 3;
    invTest->mat[1][4] = 4;
    invTest->mat[2][0] = 0;
    invTest->mat[2][1] = 0;
    invTest->mat[2][2] = 1;
    invTest->mat[2][3] = 2;
    invTest->mat[2][4] = 3;
    invTest->mat[3][0] = 0;
    invTest->mat[3][1] = 0;
    invTest->mat[3][2] = 0;
    invTest->mat[3][3] = 1;
    invTest->mat[3][4] = 2;
    invTest->mat[4][0] = 0;
    invTest->mat[4][1] = 0;
    invTest->mat[4][2] = 0;
    invTest->mat[4][3] = 0;
    invTest->mat[4][4] = 1;

    correctResult->mat[0][0] = 1;
    correctResult->mat[0][1] = -2;
    correctResult->mat[0][2] = 1;
    correctResult->mat[0][3] = 0;
    correctResult->mat[0][4] = 0;
    correctResult->mat[1][0] = 0;
    correctResult->mat[1][1] = 1;
    correctResult->mat[1][2] = -2;
    correctResult->mat[1][3] = 1;
    correctResult->mat[1][4] = 0;
    correctResult->mat[2][0] = 0;
    correctResult->mat[2][1] = 0;
    correctResult->mat[2][2] = 1;
    correctResult->mat[2][3] = -2;
    correctResult->mat[2][4] = 1;
    correctResult->mat[3][0] = 0;
    correctResult->mat[3][1] = 0;
    correctResult->mat[3][2] = 0;
    correctResult->mat[3][3] = 1;
    correctResult->mat[3][4] = -2;
    correctResult->mat[4][0] = 0;
    correctResult->mat[4][1] = 0;
    correctResult->mat[4][2] = 0;
    correctResult->mat[4][3] = 0;
    correctResult->mat[4][4] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(invTest, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(invTest);
        LOG_INFO("Inverse");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(invTest);
    FREE_MATRIX(result);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Inverse 5x5 Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseTest6x6(bool increasedLogging)
{
    struct matrix *invTest = NULL;
    struct matrix *result = NULL;
    struct matrix *correctResult = NULL;

    INIT_MATRIX(invTest, 6, 6);
    INIT_MATRIX(result, 6, 6);
    INIT_MATRIX(correctResult, 6, 6);

    LOG_INFO("Starting Inverse 6x6 Test");
    invTest->mat[0][0] = 1;
    invTest->mat[0][1] = 2;
    invTest->mat[0][2] = 3;
    invTest->mat[0][3] = 4;
    invTest->mat[0][4] = 5;
    invTest->mat[0][5] = 6;
    invTest->mat[1][0] = 0;
    invTest->mat[1][1] = 1;
    invTest->mat[1][2] = 2;
    invTest->mat[1][3] = 3;
    invTest->mat[1][4] = 4;
    invTest->mat[1][5] = 5;
    invTest->mat[2][0] = 0;
    invTest->mat[2][1] = 0;
    invTest->mat[2][2] = 1;
    invTest->mat[2][3] = 2;
    invTest->mat[2][4] = 3;
    invTest->mat[2][5] = 4;
    invTest->mat[3][0] = 0;
    invTest->mat[3][1] = 0;
    invTest->mat[3][2] = 0;
    invTest->mat[3][3] = 1;
    invTest->mat[3][4] = 2;
    invTest->mat[3][5] = 3;
    invTest->mat[4][0] = 0;
    invTest->mat[4][1] = 0;
    invTest->mat[4][2] = 0;
    invTest->mat[4][3] = 0;
    invTest->mat[4][4] = 1;
    invTest->mat[4][5] = 2;
    invTest->mat[5][0] = 0;
    invTest->mat[5][1] = 0;
    invTest->mat[5][2] = 0;
    invTest->mat[5][3] = 0;
    invTest->mat[5][4] = 0;
    invTest->mat[5][5] = 1;
    
    correctResult->mat[0][0] = 1;
    correctResult->mat[0][1] = -2;
    correctResult->mat[0][2] = 1;
    correctResult->mat[0][3] = 0;
    correctResult->mat[0][4] = 0;
    correctResult->mat[0][5] = 0;
    correctResult->mat[1][0] = 0;
    correctResult->mat[1][1] = 1;
    correctResult->mat[1][2] = -2;
    correctResult->mat[1][3] = 1;
    correctResult->mat[1][4] = 0;
    correctResult->mat[1][5] = 0;
    correctResult->mat[2][0] = 0;
    correctResult->mat[2][1] = 0;
    correctResult->mat[2][2] = 1;
    correctResult->mat[2][3] = -2;
    correctResult->mat[2][4] = 1;
    correctResult->mat[2][5] = 0;
    correctResult->mat[3][0] = 0;
    correctResult->mat[3][1] = 0;
    correctResult->mat[3][2] = 0;
    correctResult->mat[3][3] = 1;
    correctResult->mat[3][4] = -2;
    correctResult->mat[3][5] = 1;
    correctResult->mat[4][0] = 0;
    correctResult->mat[4][1] = 0;
    correctResult->mat[4][2] = 0;
    correctResult->mat[4][3] = 0;
    correctResult->mat[4][4] = 1;
    correctResult->mat[4][5] = -2;
    correctResult->mat[5][0] = 0;
    correctResult->mat[5][1] = 0;
    correctResult->mat[5][2] = 0;
    correctResult->mat[5][3] = 0;
    correctResult->mat[5][4] = 0;
    correctResult->mat[5][5] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(invTest, result));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(result, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(invTest);
        LOG_INFO("Inverse");
        printMatrix(result);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    FREE_MATRIX(invTest);
    FREE_MATRIX(result);
    FREE_MATRIX(correctResult);

    LOG_INFO("Completed Inverse 4x4 Test");

    return MATRIX_SUCCESS;
}

matrixReturnCodes copyTest(bool increasedLogging)
{
    struct matrix *copy1 = NULL;
    struct matrix *copy2 = NULL;
    struct matrix *correctResult = NULL;
    INIT_MATRIX(copy1, 2, 2);
    INIT_MATRIX(copy2, 2, 2);
    INIT_MATRIX(correctResult, 2, 2);

    //Copy test
    LOG_INFO("Starting Copy Test");
    copy1->mat[0][0] = 1;
    correctResult->mat[0][0] = 1;
    copy1->mat[0][1] = 2;
    correctResult->mat[0][1] = 2;
    copy1->mat[1][0] = 3;
    correctResult->mat[1][0] = 3;
    copy1->mat[1][1] = 4;
    correctResult->mat[1][1] = 4;

    MATRIX_MATH_RETURN_CHECK(copyMatrix(copy1, copy2));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(copy2, correctResult));

    if (increasedLogging)
    {
        LOG_INFO("Original");
        printMatrix(copy1);
        LOG_INFO("Copy");
        printMatrix(copy2);
        LOG_INFO("Correct Result");
        printMatrix(correctResult);
    }

    LOG_INFO("Completed Copy Test");

    FREE_MATRIX(copy1);
    FREE_MATRIX(copy2);

    return MATRIX_SUCCESS;
}

static matrixReturnCodes expectReturnCode(matrixReturnCodes got, matrixReturnCodes expected, const char *testName)
{
    if (got != expected)
    {
        LOG_ERROR("Expected %d but got %d in %s", expected, got, testName);
        return MATRIX_ERROR;
    }
    return MATRIX_SUCCESS;
}

static matrixReturnCodes compareMatrixApprox(struct matrix *a, struct matrix *b, matrixType eps)
{
    if (a == NULL || b == NULL)
    {
        return MATRIX_NULL_POINTER;
    }
    if (a->row != b->row || a->col != b->col)
    {
        return MATRIX_DIMENSION_MISMATCH;
    }
    for (int i = 0; i < a->row; ++i)
    {
        for (int j = 0; j < a->col; ++j)
        {
            matrixType diff = fabs(ACCESS_MATRIX(*a, i, j) - ACCESS_MATRIX(*b, i, j));
            if (diff > eps)
            {
                return MATRIX_COMPARE_FAILURE;
            }
        }
    }
    return MATRIX_SUCCESS;
}

matrixReturnCodes dimensionMismatchTests(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Dimension Mismatch Tests");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(b, 3, 3);
    INIT_MATRIX(res, 2, 2);

    MATRIX_MATH_RETURN_CHECK(expectReturnCode(addMatrix(a, b, res), MATRIX_DIMENSION_MISMATCH, "addMatrix mismatch"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(subMatrix(a, b, res), MATRIX_DIMENSION_MISMATCH, "subMatrix mismatch"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(copyMatrix(a, b), MATRIX_DIMENSION_MISMATCH, "copyMatrix mismatch"));

    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(b, 2, 2);
    INIT_MATRIX(res, 2, 2);
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(multMatrix(a, b, res), MATRIX_DIMENSION_MISMATCH, "multMatrix mismatch"));
    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(res, 2, 2);
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(transposeMatrix(a, res), MATRIX_DIMENSION_MISMATCH, "transposeMatrix mismatch"));
    FREE_MATRIX(a);
    FREE_MATRIX(res);

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(res, 2, 2);
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(scaleMatrix(a, res, 2), MATRIX_DIMENSION_MISMATCH, "scaleMatrix mismatch"));
    FREE_MATRIX(a);
    FREE_MATRIX(res);

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(res, 2, 3);
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(inverseMatrix(a, res), MATRIX_DIMENSION_MISMATCH, "inverseMatrix non-square"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(identityMatrixMinusA(a, res), MATRIX_DIMENSION_MISMATCH, "identityMatrixMinusA non-square"));
    FREE_MATRIX(a);
    FREE_MATRIX(res);

    LOG_INFO("Completed Dimension Mismatch Tests");
    return MATRIX_SUCCESS;
}

matrixReturnCodes nanAndInfTests(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting NaN/Inf Tests");

    struct matrix *a = NULL;
    INIT_MATRIX(a, 2, 2);

    a->mat[0][0] = NAN;
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(nanCheckMatrix(a), MATRIX_NAN_FAILURE, "nanCheckMatrix NaN"));

    a->mat[0][0] = INFINITY;
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(nanCheckMatrix(a), MATRIX_INF_FAILURE, "nanCheckMatrix Inf"));

    FREE_MATRIX(a);

    LOG_INFO("Completed NaN/Inf Tests");
    return MATRIX_SUCCESS;
}

matrixReturnCodes staticMatrixAddTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Static Matrix Add Test");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    STATIC_MATRIX_DIRECTIVE(a, 2, 2, a);
    STATIC_MATRIX_DIRECTIVE(b, 2, 2, b);
    STATIC_MATRIX_DIRECTIVE(res, 2, 2, res);
    STATIC_MATRIX_DIRECTIVE(correct, 2, 2, correct);

    ACCESS_STATIC_MATRIX(*a, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*a, 0, 1) = 2;
    ACCESS_STATIC_MATRIX(*a, 1, 0) = 3;
    ACCESS_STATIC_MATRIX(*a, 1, 1) = 4;

    ACCESS_STATIC_MATRIX(*b, 0, 0) = 5;
    ACCESS_STATIC_MATRIX(*b, 0, 1) = 6;
    ACCESS_STATIC_MATRIX(*b, 1, 0) = 7;
    ACCESS_STATIC_MATRIX(*b, 1, 1) = 8;

    ACCESS_STATIC_MATRIX(*correct, 0, 0) = 6;
    ACCESS_STATIC_MATRIX(*correct, 0, 1) = 8;
    ACCESS_STATIC_MATRIX(*correct, 1, 0) = 10;
    ACCESS_STATIC_MATRIX(*correct, 1, 1) = 12;

    MATRIX_MATH_RETURN_CHECK(addMatrix(a, b, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    LOG_INFO("Completed Static Matrix Add Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes staticMatrixMultTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Static Matrix Mult Test");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    STATIC_MATRIX_DIRECTIVE(a, 2, 2, a);
    STATIC_MATRIX_DIRECTIVE(b, 2, 2, b);
    STATIC_MATRIX_DIRECTIVE(res, 2, 2, res);
    STATIC_MATRIX_DIRECTIVE(correct, 2, 2, correct);

    ACCESS_STATIC_MATRIX(*a, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*a, 0, 1) = 2;
    ACCESS_STATIC_MATRIX(*a, 1, 0) = 3;
    ACCESS_STATIC_MATRIX(*a, 1, 1) = 4;

    ACCESS_STATIC_MATRIX(*b, 0, 0) = 5;
    ACCESS_STATIC_MATRIX(*b, 0, 1) = 6;
    ACCESS_STATIC_MATRIX(*b, 1, 0) = 7;
    ACCESS_STATIC_MATRIX(*b, 1, 1) = 8;

    ACCESS_STATIC_MATRIX(*correct, 0, 0) = 19;
    ACCESS_STATIC_MATRIX(*correct, 0, 1) = 22;
    ACCESS_STATIC_MATRIX(*correct, 1, 0) = 43;
    ACCESS_STATIC_MATRIX(*correct, 1, 1) = 50;

    MATRIX_MATH_RETURN_CHECK(multMatrix(a, b, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    LOG_INFO("Completed Static Matrix Mult Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes staticMatrixTransposeTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Static Matrix Transpose Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    STATIC_MATRIX_DIRECTIVE(a, 2, 3, a);
    STATIC_MATRIX_DIRECTIVE(res, 3, 2, res);
    STATIC_MATRIX_DIRECTIVE(correct, 3, 2, correct);

    ACCESS_STATIC_MATRIX(*a, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*a, 0, 1) = 2;
    ACCESS_STATIC_MATRIX(*a, 0, 2) = 3;
    ACCESS_STATIC_MATRIX(*a, 1, 0) = 4;
    ACCESS_STATIC_MATRIX(*a, 1, 1) = 5;
    ACCESS_STATIC_MATRIX(*a, 1, 2) = 6;

    ACCESS_STATIC_MATRIX(*correct, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*correct, 0, 1) = 4;
    ACCESS_STATIC_MATRIX(*correct, 1, 0) = 2;
    ACCESS_STATIC_MATRIX(*correct, 1, 1) = 5;
    ACCESS_STATIC_MATRIX(*correct, 2, 0) = 3;
    ACCESS_STATIC_MATRIX(*correct, 2, 1) = 6;

    MATRIX_MATH_RETURN_CHECK(transposeMatrix(a, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    LOG_INFO("Completed Static Matrix Transpose Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes staticMatrixScaleTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Static Matrix Scale Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    STATIC_MATRIX_DIRECTIVE(a, 2, 2, a);
    STATIC_MATRIX_DIRECTIVE(res, 2, 2, res);
    STATIC_MATRIX_DIRECTIVE(correct, 2, 2, correct);

    ACCESS_STATIC_MATRIX(*a, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*a, 0, 1) = 2;
    ACCESS_STATIC_MATRIX(*a, 1, 0) = 3;
    ACCESS_STATIC_MATRIX(*a, 1, 1) = 4;

    ACCESS_STATIC_MATRIX(*correct, 0, 0) = 3;
    ACCESS_STATIC_MATRIX(*correct, 0, 1) = 6;
    ACCESS_STATIC_MATRIX(*correct, 1, 0) = 9;
    ACCESS_STATIC_MATRIX(*correct, 1, 1) = 12;

    MATRIX_MATH_RETURN_CHECK(scaleMatrix(a, res, 3));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    LOG_INFO("Completed Static Matrix Scale Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes staticMatrixIdentityMinusATest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Static Matrix I - A Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    STATIC_MATRIX_DIRECTIVE(a, 2, 2, a);
    STATIC_MATRIX_DIRECTIVE(res, 2, 2, res);
    STATIC_MATRIX_DIRECTIVE(correct, 2, 2, correct);

    ACCESS_STATIC_MATRIX(*a, 0, 0) = 1;
    ACCESS_STATIC_MATRIX(*a, 0, 1) = 2;
    ACCESS_STATIC_MATRIX(*a, 1, 0) = 3;
    ACCESS_STATIC_MATRIX(*a, 1, 1) = 4;

    ACCESS_STATIC_MATRIX(*correct, 0, 0) = 0;
    ACCESS_STATIC_MATRIX(*correct, 0, 1) = -2;
    ACCESS_STATIC_MATRIX(*correct, 1, 0) = -3;
    ACCESS_STATIC_MATRIX(*correct, 1, 1) = -3;

    MATRIX_MATH_RETURN_CHECK(identityMatrixMinusA(a, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    LOG_INFO("Completed Static Matrix I - A Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes compareDimensionMismatchTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Compare Dimension Mismatch Test");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(b, 2, 3);

    MATRIX_MATH_RETURN_CHECK(expectReturnCode(compareMatrieces(a, b), MATRIX_DIMENSION_MISMATCH, "compareMatrieces mismatch"));

    FREE_MATRIX(a);
    FREE_MATRIX(b);

    LOG_INFO("Completed Compare Dimension Mismatch Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes copy2DArrayTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting 2D Array Copy Test");

    ekfType data[2][3] = {{1, 2, 3}, {4, 5, 6}};

    struct matrix *m = NULL;
    INIT_MATRIX(m, 2, 3);
    COPY_2DARRAY_TO_MATRIX(data, m);

    if (ACCESS_MATRIX(*m, 0, 0) != 1 || ACCESS_MATRIX(*m, 0, 1) != 2 || ACCESS_MATRIX(*m, 0, 2) != 3 ||
        ACCESS_MATRIX(*m, 1, 0) != 4 || ACCESS_MATRIX(*m, 1, 1) != 5 || ACCESS_MATRIX(*m, 1, 2) != 6)
    {
        LOG_ERROR("COPY_2DARRAY_TO_MATRIX failed to copy values to malloc matrix.");
        FREE_MATRIX(m);
        return MATRIX_ERROR;
    }

    struct matrix *s = NULL;
    STATIC_MATRIX_DIRECTIVE(s, 2, 3, s);
    COPY_2DARRAY_TO_MATRIX(data, s);

    if (ACCESS_MATRIX(*s, 0, 0) != 1 || ACCESS_MATRIX(*s, 0, 1) != 2 || ACCESS_MATRIX(*s, 0, 2) != 3 ||
        ACCESS_MATRIX(*s, 1, 0) != 4 || ACCESS_MATRIX(*s, 1, 1) != 5 || ACCESS_MATRIX(*s, 1, 2) != 6)
    {
        LOG_ERROR("COPY_2DARRAY_TO_MATRIX failed to copy values to static matrix.");
        FREE_MATRIX(m);
        return MATRIX_ERROR;
    }

    FREE_MATRIX(m);
    LOG_INFO("Completed 2D Array Copy Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes multMatrixZeroingTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Mult Matrix Zeroing Test");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(b, 2, 2);
    INIT_MATRIX(res, 2, 2);

    a->mat[0][0] = 0;
    a->mat[0][1] = 0;
    a->mat[1][0] = 0;
    a->mat[1][1] = 0;

    b->mat[0][0] = 1;
    b->mat[0][1] = 2;
    b->mat[1][0] = 3;
    b->mat[1][1] = 4;

    res->mat[0][0] = 9;
    res->mat[0][1] = 9;
    res->mat[1][0] = 9;
    res->mat[1][1] = 9;

    MATRIX_MATH_RETURN_CHECK(multMatrix(a, b, res));

    if (res->mat[0][0] != 0 || res->mat[0][1] != 0 || res->mat[1][0] != 0 || res->mat[1][1] != 0)
    {
        LOG_ERROR("multMatrix did not zero the result matrix before multiplying.");
        FREE_MATRIX(a);
        FREE_MATRIX(b);
        FREE_MATRIX(res);
        return MATRIX_ERROR;
    }

    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);

    LOG_INFO("Completed Mult Matrix Zeroing Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseSingularTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Singular Inverse Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(res, 2, 2);

    a->mat[0][0] = 1;
    a->mat[0][1] = 2;
    a->mat[1][0] = 2;
    a->mat[1][1] = 4;

    matrixReturnCodes rc = inverseMatrix(a, res);
    if (rc != MATRIX_NAN_FAILURE && rc != MATRIX_INF_FAILURE)
    {
        LOG_ERROR("inverseMatrix did not flag singular matrix with NAN/INF return code.");
        FREE_MATRIX(a);
        FREE_MATRIX(res);
        return MATRIX_ERROR;
    }

    FREE_MATRIX(a);
    FREE_MATRIX(res);

    LOG_INFO("Completed Singular Inverse Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseJitterTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Inverse Jitter Test");

    struct matrix *a = NULL;
    struct matrix *inv = NULL;
    struct matrix *res = NULL;
    struct matrix *identity = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(inv, 2, 2);
    INIT_MATRIX(res, 2, 2);
    INIT_MATRIX(identity, 2, 2);

    a->mat[0][0] = 1;
    a->mat[0][1] = 1;
    a->mat[1][0] = 1;
    a->mat[1][1] = 1;

    identity->mat[0][0] = 1;
    identity->mat[0][1] = 0;
    identity->mat[1][0] = 0;
    identity->mat[1][1] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrixWithJitter(a, inv, (matrixType)1e-6, 3, (matrixType)100));
    MATRIX_MATH_RETURN_CHECK(multMatrix(a, inv, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrixApprox(res, identity, 1e-3));

    if (a->mat[0][0] < (matrixType)(1 + 1e-6) || a->mat[1][1] < (matrixType)(1 + 1e-6))
    {
        LOG_ERROR("inverseMatrixWithJitter did not add jitter to the diagonal.");
        FREE_MATRIX(a);
        FREE_MATRIX(inv);
        FREE_MATRIX(res);
        FREE_MATRIX(identity);
        return MATRIX_ERROR;
    }

    FREE_MATRIX(a);
    FREE_MATRIX(inv);
    FREE_MATRIX(res);
    FREE_MATRIX(identity);

    LOG_INFO("Completed Inverse Jitter Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes nanPropagationTests(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting NaN Propagation Tests");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(b, 2, 2);
    INIT_MATRIX(res, 2, 2);

    a->mat[0][0] = NAN;
    b->mat[0][0] = 1;

    MATRIX_MATH_RETURN_CHECK(expectReturnCode(addMatrix(a, b, res), MATRIX_NAN_FAILURE, "addMatrix NaN input"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(subMatrix(a, b, res), MATRIX_NAN_FAILURE, "subMatrix NaN input"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(multMatrix(a, b, res), MATRIX_NAN_FAILURE, "multMatrix NaN input"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(scaleMatrix(a, res, 2), MATRIX_NAN_FAILURE, "scaleMatrix NaN input"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(transposeMatrix(a, res), MATRIX_NAN_FAILURE, "transposeMatrix NaN input"));

    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);

    LOG_INFO("Completed NaN Propagation Tests");
    return MATRIX_SUCCESS;
}

matrixReturnCodes identityPropertiesTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Identity Properties Test");

    struct matrix *a = NULL;
    struct matrix *i = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    INIT_MATRIX(a, 3, 3);
    INIT_MATRIX(i, 3, 3);
    INIT_MATRIX(res, 3, 3);
    INIT_MATRIX(correct, 3, 3);

    a->mat[0][0] = 1; a->mat[0][1] = 2; a->mat[0][2] = 3;
    a->mat[1][0] = 4; a->mat[1][1] = 5; a->mat[1][2] = 6;
    a->mat[2][0] = 7; a->mat[2][1] = 8; a->mat[2][2] = 9;

    for (int irow = 0; irow < 3; ++irow)
    {
        for (int jcol = 0; jcol < 3; ++jcol)
        {
            i->mat[irow][jcol] = (irow == jcol) ? 1 : 0;
        }
    }

    MATRIX_MATH_RETURN_CHECK(multMatrix(a, i, res));
    MATRIX_MATH_RETURN_CHECK(copyMatrix(a, correct));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    MATRIX_MATH_RETURN_CHECK(multMatrix(i, a, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    FREE_MATRIX(a);
    FREE_MATRIX(i);
    FREE_MATRIX(res);
    FREE_MATRIX(correct);

    LOG_INFO("Completed Identity Properties Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes setIdentityMatrixTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Set Identity Matrix Test");

    struct matrix *a = NULL;
    struct matrix *expected = NULL;

    INIT_MATRIX(a, 3, 3);
    INIT_MATRIX(expected, 3, 3);

    a->mat[0][0] = 1; a->mat[0][1] = 2; a->mat[0][2] = 3;
    a->mat[1][0] = 4; a->mat[1][1] = 5; a->mat[1][2] = 6;
    a->mat[2][0] = 7; a->mat[2][1] = 8; a->mat[2][2] = 9;

    expected->mat[0][0] = 1; expected->mat[0][1] = 0; expected->mat[0][2] = 0;
    expected->mat[1][0] = 0; expected->mat[1][1] = 1; expected->mat[1][2] = 0;
    expected->mat[2][0] = 0; expected->mat[2][1] = 0; expected->mat[2][2] = 1;

    MATRIX_MATH_RETURN_CHECK(setIdentityMatrix(a));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(a, expected));

    FREE_MATRIX(a);
    FREE_MATRIX(expected);

    LOG_INFO("Completed Set Identity Matrix Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes transposeTwiceTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Transpose Twice Test");

    struct matrix *a = NULL;
    struct matrix *t1 = NULL;
    struct matrix *t2 = NULL;

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(t1, 3, 2);
    INIT_MATRIX(t2, 2, 3);

    a->mat[0][0] = 1; a->mat[0][1] = 2; a->mat[0][2] = 3;
    a->mat[1][0] = 4; a->mat[1][1] = 5; a->mat[1][2] = 6;

    MATRIX_MATH_RETURN_CHECK(transposeMatrix(a, t1));
    MATRIX_MATH_RETURN_CHECK(transposeMatrix(t1, t2));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(a, t2));

    FREE_MATRIX(a);
    FREE_MATRIX(t1);
    FREE_MATRIX(t2);

    LOG_INFO("Completed Transpose Twice Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes scaleIdentityZeroTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Scale Identity/Zero Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(res, 2, 2);
    INIT_MATRIX(correct, 2, 2);

    a->mat[0][0] = 1; a->mat[0][1] = -2;
    a->mat[1][0] = 3; a->mat[1][1] = 4;

    MATRIX_MATH_RETURN_CHECK(scaleMatrix(a, res, 1));
    MATRIX_MATH_RETURN_CHECK(copyMatrix(a, correct));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    MATRIX_MATH_RETURN_CHECK(scaleMatrix(a, res, 0));
    if (res->mat[0][0] != 0 || res->mat[0][1] != 0 || res->mat[1][0] != 0 || res->mat[1][1] != 0)
    {
        LOG_ERROR("scaleMatrix with 0 did not produce all zeros.");
        FREE_MATRIX(a);
        FREE_MATRIX(res);
        FREE_MATRIX(correct);
        return MATRIX_ERROR;
    }

    FREE_MATRIX(a);
    FREE_MATRIX(res);
    FREE_MATRIX(correct);

    LOG_INFO("Completed Scale Identity/Zero Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseSanityTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Inverse Sanity Test");

    struct matrix *a = NULL;
    struct matrix *inv = NULL;
    struct matrix *res = NULL;
    struct matrix *identity = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(inv, 2, 2);
    INIT_MATRIX(res, 2, 2);
    INIT_MATRIX(identity, 2, 2);

    a->mat[0][0] = 4; a->mat[0][1] = 7;
    a->mat[1][0] = 2; a->mat[1][1] = 6;

    identity->mat[0][0] = 1; identity->mat[0][1] = 0;
    identity->mat[1][0] = 0; identity->mat[1][1] = 1;

    MATRIX_MATH_RETURN_CHECK(inverseMatrix(a, inv));
    MATRIX_MATH_RETURN_CHECK(multMatrix(a, inv, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrixApprox(res, identity, 1e-6));

    FREE_MATRIX(a);
    FREE_MATRIX(inv);
    FREE_MATRIX(res);
    FREE_MATRIX(identity);

    LOG_INFO("Completed Inverse Sanity Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes inverseDoesNotMutateTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Inverse Non-Mutate Test");

    struct matrix *a = NULL;
    struct matrix *orig = NULL;
    struct matrix *inv = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(orig, 2, 2);
    INIT_MATRIX(inv, 2, 2);

    a->mat[0][0] = 1; a->mat[0][1] = 2;
    a->mat[1][0] = 3; a->mat[1][1] = 4;

    MATRIX_MATH_RETURN_CHECK(copyMatrix(a, orig));
    MATRIX_MATH_RETURN_CHECK(inverseMatrix(a, inv));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(a, orig));

    FREE_MATRIX(a);
    FREE_MATRIX(orig);
    FREE_MATRIX(inv);

    LOG_INFO("Completed Inverse Non-Mutate Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes jaggedStaticParityTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Jagged vs Static Parity Test");

    struct matrix *j = NULL;
    struct matrix *jres = NULL;
    struct matrix *s = NULL;
    struct matrix *sres = NULL;

    INIT_MATRIX(j, 2, 3);
    INIT_MATRIX(jres, 3, 2);
    STATIC_MATRIX_DIRECTIVE(s, 2, 3, s);
    STATIC_MATRIX_DIRECTIVE(sres, 3, 2, sres);

    j->mat[0][0] = 1; j->mat[0][1] = 2; j->mat[0][2] = 3;
    j->mat[1][0] = 4; j->mat[1][1] = 5; j->mat[1][2] = 6;

    ekfType data[2][3] = {{1, 2, 3}, {4, 5, 6}};
    COPY_2DARRAY_TO_MATRIX(data, s);

    MATRIX_MATH_RETURN_CHECK(transposeMatrix(j, jres));
    MATRIX_MATH_RETURN_CHECK(transposeMatrix(s, sres));

    MATRIX_MATH_RETURN_CHECK(compareMatrixApprox(jres, sres, 0));

    FREE_MATRIX(j);
    FREE_MATRIX(jres);

    LOG_INFO("Completed Jagged vs Static Parity Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes nonsquareMultiplyTest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Non-square Multiply Test");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;
    struct matrix *correct = NULL;

    INIT_MATRIX(a, 2, 3);
    INIT_MATRIX(b, 3, 2);
    INIT_MATRIX(res, 2, 2);
    INIT_MATRIX(correct, 2, 2);

    a->mat[0][0] = 1; a->mat[0][1] = 2; a->mat[0][2] = 3;
    a->mat[1][0] = 4; a->mat[1][1] = 5; a->mat[1][2] = 6;

    b->mat[0][0] = 7;  b->mat[0][1] = 8;
    b->mat[1][0] = 9;  b->mat[1][1] = 10;
    b->mat[2][0] = 11; b->mat[2][1] = 12;

    correct->mat[0][0] = 58;
    correct->mat[0][1] = 64;
    correct->mat[1][0] = 139;
    correct->mat[1][1] = 154;

    MATRIX_MATH_RETURN_CHECK(multMatrix(a, b, res));
    MATRIX_MATH_RETURN_CHECK(compareMatrieces(res, correct));

    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);
    FREE_MATRIX(correct);

    LOG_INFO("Completed Non-square Multiply Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes nanIdentityMinusATest(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting NaN IdentityMinusA Test");

    struct matrix *a = NULL;
    struct matrix *res = NULL;
    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(res, 2, 2);

    a->mat[0][0] = NAN;
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(identityMatrixMinusA(a, res), MATRIX_NAN_FAILURE, "identityMatrixMinusA NaN input"));

    FREE_MATRIX(a);
    FREE_MATRIX(res);

    LOG_INFO("Completed NaN IdentityMinusA Test");
    return MATRIX_SUCCESS;
}

matrixReturnCodes aliasingGuardTests(bool increasedLogging)
{
    (void)increasedLogging;
    LOG_INFO("Starting Aliasing Guard Tests");

    struct matrix *a = NULL;
    struct matrix *b = NULL;
    struct matrix *res = NULL;

    INIT_MATRIX(a, 2, 2);
    INIT_MATRIX(b, 2, 2);
    INIT_MATRIX(res, 2, 2);

    // multMatrix: res must not alias a or b
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(multMatrix(a, b, a), MATRIX_ERROR, "multMatrix alias a"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(multMatrix(a, b, b), MATRIX_ERROR, "multMatrix alias b"));

    // addMatrix: res must not alias a or b
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(addMatrix(a, b, a), MATRIX_ERROR, "addMatrix alias a"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(addMatrix(a, b, b), MATRIX_ERROR, "addMatrix alias b"));

    // subMatrix: res must not alias a or b
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(subMatrix(a, b, a), MATRIX_ERROR, "subMatrix alias a"));
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(subMatrix(a, b, b), MATRIX_ERROR, "subMatrix alias b"));

    // transposeMatrix: res must not alias a
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(transposeMatrix(a, a), MATRIX_ERROR, "transposeMatrix alias a"));

    // inverseMatrix: res must not alias a
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(inverseMatrix(a, a), MATRIX_ERROR, "inverseMatrix alias a"));

    // inverseMatrixWithJitter: res must not alias a
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(inverseMatrixWithJitter(a, a, (matrixType)1e-6, 3, (matrixType)100), MATRIX_ERROR, "inverseMatrixWithJitter alias a"));

    // identityMatrixMinusA: res must not alias a
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(identityMatrixMinusA(a, a), MATRIX_ERROR, "identityMatrixMinusA alias a"));

    // copyMatrix: res must not alias a
    MATRIX_MATH_RETURN_CHECK(expectReturnCode(copyMatrix(a, a), MATRIX_ERROR, "copyMatrix alias a"));

    // scaleMatrix: in-place is allowed
    MATRIX_MATH_RETURN_CHECK(scaleMatrix(a, a, 2));

    FREE_MATRIX(a);
    FREE_MATRIX(b);
    FREE_MATRIX(res);

    LOG_INFO("Completed Aliasing Guard Tests");
    return MATRIX_SUCCESS;
}
