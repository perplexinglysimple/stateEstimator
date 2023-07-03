/**
 * @file ekfMatrixMathTest.c
 * 
 * I know this may leak memory if malloc fails and all other things... Its a test. I might fix it I might not.
*/

#include "EKF.h"

matrixReturnCodes copyTest(bool increasedLogging);
matrixReturnCodes inverseTest(bool increasedLogging);
matrixReturnCodes additionTest(bool increasedLogging);
matrixReturnCodes subtractionTest(bool increasedLogging);
matrixReturnCodes multiplicationTest(bool increasedLogging);
matrixReturnCodes tranposeTest(bool increasedLogging);
matrixReturnCodes scaleTest(bool increasedLogging);
matrixReturnCodes idenityMatrixMinusATest(bool increasedLogging);

int main()
{
    LOG_INFO("Starting Matrix Math Tests");
    bool increasedLogging = true;

    // Copy Test
    MATRIX_MATH_RETURN_CHECK(copyTest(increasedLogging));

    //Inverse test
    MATRIX_MATH_RETURN_CHECK(inverseTest(increasedLogging));

    //Add test
    MATRIX_MATH_RETURN_CHECK(additionTest(increasedLogging));

    //Sub test
    MATRIX_MATH_RETURN_CHECK(subtractionTest(increasedLogging));

    //Mult test
    MATRIX_MATH_RETURN_CHECK(multiplicationTest(increasedLogging));

    //Transpose test
    MATRIX_MATH_RETURN_CHECK(tranposeTest(increasedLogging));

    //Scale test
    MATRIX_MATH_RETURN_CHECK(scaleTest(increasedLogging));

    //I - A test
    MATRIX_MATH_RETURN_CHECK(idenityMatrixMinusATest(increasedLogging));

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