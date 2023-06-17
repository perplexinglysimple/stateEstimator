//#define TEST_MULT

#include "matrixMath.h"

//Inputs are three matrixes. a and b matrix will be multiplied and the result will be output to res
matrixReturnCodes multMatrix(struct matrix *a, struct matrix *b, struct matrix *res) {
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX(b);
	
	NULL_CHECK_MATRIX_RES(res);
	
	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(b);
	NON_INIT_CHECK_MATRIX(res);
	
	DIMENSION_CHECK_MULT_MATRIX(a, b, res);
	
	int arow, acol, bcol;
	for(arow = 0; arow < a->row; ++arow) {
		for(bcol = 0; bcol < b->col; ++bcol) {
			for(acol = 0; acol < a->col; ++acol) {
				res->mat[arow][bcol] += a->mat[arow][acol] * b->mat[acol][bcol];
				#ifdef VERBOSE_TEST
					printMatrix(res);
					printf("  arow=%d, acol=%d, bcol=%d\n", arow, acol, bcol);
					printf("  aval=%d, bval=%d\n", a->mat[arow][acol], b->mat[acol][bcol]);
				#endif //VERBOSE_TEST
			}
		}
	}
	return MATRIX_SUCCESS;
}

matrixReturnCodes scaleMatrix(struct matrix *a, struct matrix *res, matrixType scaler) {
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX_RES(res);

	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(res);
	
	DIMENSION_CHECK_SCALER_MATRIX(a, res);
	
	int i, j;
	for(i = 0; i < a->row; ++i) {
		for(j = 0; j < a->col; ++j) {
			res->mat[i][j] = a->mat[i][j] * scaler;
		}
	}
	return MATRIX_SUCCESS;
}

matrixReturnCodes addMatrix(struct matrix *a, struct matrix *b, struct matrix *res) {
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX(b);
	
	NULL_CHECK_MATRIX_RES(res);
	
	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(b);
	NON_INIT_CHECK_MATRIX(res);
	
	DIMENSION_CHECK_ADD_MATRIX(a, b, res);

	int arow, acol;
	for(arow = 0; arow < a->row; ++arow) {
		for(acol = 0; acol < a->col; ++acol) {
			res->mat[arow][acol] = a->mat[arow][acol] + b->mat[arow][acol];
		}
	}
	
	return MATRIX_SUCCESS;
}

void printMatrix(struct matrix *a) {
	if(a == NULL || !a->initilized) {
		printf("Trying to print a uninitilized matrix");
		return;
	}
	int i, j;
	for(i = 0; i < a->row; ++i) {
		if(i > 0)
			printf("\n");
		for(j = 0; j < a->col; ++j) {
			#ifdef _INT
				printf("%d, ", a->mat[i][j]);
			#elif defined(_DOUBLE)
				printf("%f, ", a->mat[i][j]);
			#else
				printf("           Type not defined!");
				printf("%d, ", a->mat[i][j]);
			#endif
		}
	}
	printf("\n");
}

matrixReturnCodes compareMatrieces(struct matrix *a, struct matrix *b) {
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX(b);

	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(b);
	if(a->row != b->row || a->col != a->col) {
		return MATRIX_DIMENSION_MISMATCH;
	}
	int i, j;
	for(i = 0; i < a->row; ++i) {
		for(j = 0; j < a->col; ++j) {
			if(a->mat[i][j] != b->mat[i][j]) {
				return MATRIX_COMPARE_FAILURE;
			}
		}
	}
	return MATRIX_SUCCESS;
}

matrixReturnCodes tranposeMatrix(struct matrix *a, struct matrix *b) {
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX(b);

	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(b);
	if(a->row != b->col || a->col != b->row) {
		return MATRIX_DIMENSION_MISMATCH;
	}
	int i, j;
	for(i = 0; i < a->row; ++i) {
		for(j = 0; j < a->col; ++j) {
			b->mat[j][i] = a->mat[i][j];
		}
	}
	return MATRIX_SUCCESS;
}

matrixReturnCodes subMatrix(struct matrix *a, struct matrix *b, struct matrix *res)
{
	NULL_CHECK_MATRIX(a);
	NULL_CHECK_MATRIX(b);
	
	NULL_CHECK_MATRIX_RES(res);
	
	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(b);
	NON_INIT_CHECK_MATRIX(res);
	
	DIMENSION_CHECK_ADD_MATRIX(a, b, res);

	int arow, acol;
	for(arow = 0; arow < a->row; ++arow) {
		for(acol = 0; acol < a->col; ++acol) {
			res->mat[arow][acol] = a->mat[arow][acol] - b->mat[arow][acol];
		}
	}
	
	return MATRIX_SUCCESS;
}

matrixReturnCodes inverseMatrix(struct matrix *a, struct matrix *res)
{
	NULL_CHECK_MATRIX(a);
	
	NULL_CHECK_MATRIX_RES(res);
	
	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(res);

	LOG_ERROR("inverseMatrix not implemented");
	return MATRIX_ERROR;
}

matrixReturnCodes idenityMatrixMinusA(struct matrix *a, struct matrix *res)
{
	return MATRIX_ERROR;
}

matrixReturnCodes copyMatrix(struct matrix *a, struct matrix *res)
{
	NULL_CHECK_MATRIX(a);

	NULL_CHECK_MATRIX_RES(res);

	NON_INIT_CHECK_MATRIX(a);
	NON_INIT_CHECK_MATRIX(res);

	if(a->row != res->row || a->col != res->col) {
		return MATRIX_DIMENSION_MISMATCH;
	}
	int i, j;
	for(i = 0; i < a->row; ++i) {
		for(j = 0; j < a->col; ++j) {
			res->mat[i][j] = a->mat[i][j];
		}
	}
	return MATRIX_SUCCESS;
}

//Test for all this shit
#ifdef TEST_MULT

int main() {
	struct matrix *A = NULL;
	struct matrix *B = NULL;
	struct matrix *RES = NULL;
	struct matrix *REALRES = NULL;
	printf("Here \n");
	INIT_MATRIX(A, 3, 3);
	printf("After A INIT\n");
	INIT_MATRIX(B, 3, 2);
	printf("After B INIT\n");
	INIT_MATRIX(RES, 3, 2);
	printf("After RES INIT\n");
	
	A->mat[0][0] = 6;
	A->mat[0][1] = 3;
	A->mat[0][2] = 0;
	
	A->mat[1][0] = 2;
	A->mat[1][1] = 5;
	A->mat[1][2] = 1;
	
	A->mat[2][0] = 9;
	A->mat[2][1] = 8;
	A->mat[2][2] = 6;
	
	printf("Printing A \n");
	printMatrix(A);
	printf("\n");
	
	B->mat[0][0] = 7;
	B->mat[0][1] = 4;
	
	B->mat[1][0] = 6;
	B->mat[1][1] = 7;
	
	B->mat[2][0] = 5;
	B->mat[2][1] = 0;
	
	printf("Printing B \n");
	printMatrix(B);
	printf("\n");
	
	int returnVal = multMatrix(A, B, RES);
	
	if(returnVal != 0) {
		//print something
		printf("Failed with a code %d\n", returnVal);
	}
	else {
		//print the resulting matrix
		printf("Printing multiplication result \n");
		printMatrix(RES);
		printf("\n");
		INIT_MATRIX(REALRES, 3, 2);
		REALRES->mat[0][0] = 60;
	    REALRES->mat[0][1] = 45;
		REALRES->mat[1][0] = 49;
		REALRES->mat[1][1] = 43;
		REALRES->mat[2][0] = 141;
		REALRES->mat[2][1] = 92;
		printf("Printing right multiplication result \n");
		printMatrix(REALRES);
		printf("\n");
		printf("Comparing results... and it %s \n", compareMatrieces(REALRES,RES)==0?"matches":"does not match");
	}
	assert(returnVal == 0);
	FREE_MATRIX(A);
	FREE_MATRIX(B);
	FREE_MATRIX(RES);
	FREE_MATRIX(REALRES);
	
	printf("Printing Test Two \n");
	printf("Here \n");
	INIT_MATRIX(A, 4, 4);
	printf("After A \n");
	INIT_MATRIX(B, 4, 1);
	printf("After B \n");
	INIT_MATRIX(RES, 4, 1);
	printf("After RES \n");
	
	A->mat[0][0] = 1;
	//A->mat[0][1] = 0;
	A->mat[0][2] = 2;
	//A->mat[0][3] = 0;
	
	//A->mat[1][0] = 0;
	A->mat[1][1] = 3;
	//A->mat[1][2] = 0;
	A->mat[1][3] = 4;
	
	//A->mat[2][0] = 0;
	//A->mat[2][1] = 0;
	A->mat[2][2] = 5;
	//A->mat[2][3] = 0;
	
	A->mat[3][0] = 6;
	//A->mat[3][1] = 0;
	//A->mat[3][2] = 0;
	A->mat[3][3] = 7;
	
	printf("Printing A \n");
	printMatrix(A);
	printf("\n");
	
	B->mat[0][0] = 2;
	B->mat[1][0] = 5;
	B->mat[2][0] = 1;
	B->mat[3][0] = 8;
	printf("Printing B \n");
	printMatrix(B);
	printf("\n");
	
	returnVal = multMatrix(A, B, RES);
	
	if(returnVal != 0) {
		//print something
		printf("Failed with a code %d\n", returnVal);
	}
	else {
		//print the resulting matrix
		printf("Printing multiplication result \n");
		printMatrix(RES);
		printf("\n");
		INIT_MATRIX(REALRES, 4, 1);
		REALRES->mat[0][0] = 4;
		REALRES->mat[1][0] = 47;
		REALRES->mat[2][0] = 5;
		REALRES->mat[3][0] = 68;
		printf("Printing right multiplication result \n");
		printMatrix(REALRES);
		printf("\n");
		printf("Comparing results... and it %s \n", compareMatrieces(REALRES,RES)==0?"matches":"does not match");
	}
	assert(returnVal == 0);
	FREE_MATRIX(A);
	FREE_MATRIX(B);
	FREE_MATRIX(RES);
	FREE_MATRIX(REALRES);
	return 0;
}

#endif //TEST_MULT

#ifdef TEST_ADD
int main() {
	
	struct matrix *A = NULL;
	struct matrix *B = NULL;
	struct matrix *RES = NULL;
	struct matrix *REALRES = NULL;
	printf("Here \n");
	INIT_MATRIX(A, 3, 3);
	printf("After A \n");
	INIT_MATRIX(B, 3, 3);
	printf("After B \n");
	INIT_MATRIX(RES, 3, 3);
	printf("After RES \n");

	A->mat[0][0] = 6;
	A->mat[0][1] = 3;
	A->mat[0][2] = 0;
	
	A->mat[1][0] = 2;
	A->mat[1][1] = 5;
	A->mat[1][2] = 1;
	
	A->mat[2][0] = 9;
	A->mat[2][1] = 8;
	A->mat[2][2] = 6;
	
	printf("Printing A \n");
	printMatrix(A);
	printf("\n");
	
	B->mat[0][0] = 7;
	B->mat[0][1] = 4;
	B->mat[0][2] = 4;
	
	B->mat[1][0] = 6;
	B->mat[1][1] = 7;
	B->mat[1][2] = 4;
	
	B->mat[2][0] = 5;
	B->mat[2][1] = 0;
	B->mat[2][2] = 4;
	
	printf("Printing B \n");
	printMatrix(B);
	printf("\n");
	
	int returnVal = addMatrix(A, B, RES);
	
	if(returnVal != 0) {
		//print something
		printf("Failed with a code %d\n", returnVal);
	}
	else {
		//print the resulting matrix
		printf("Printing addition result \n");
		printMatrix(RES);
		printf("\n");
		INIT_MATRIX(REALRES, 3, 3);
		REALRES->mat[0][0] = 13;
		REALRES->mat[0][1] = 7;
		REALRES->mat[0][2] = 4;
		
		REALRES->mat[1][0] = 8;
		REALRES->mat[1][1] = 12;
		REALRES->mat[1][2] = 5;
		
		REALRES->mat[2][0] = 14;
		REALRES->mat[2][1] = 8;
		REALRES->mat[2][2] = 10;
		printf("Printing right addition result \n");
		printMatrix(REALRES);
		printf("\n");
		printf("Comparing results... and it %s \n", compareMatrieces(REALRES,RES)==0?"matches":"does not match");
	}
	assert(returnVal == 0);
	FREE_MATRIX(A);
	FREE_MATRIX(B);
	FREE_MATRIX(RES);
	FREE_MATRIX(REALRES);
	return 0;
}
#endif //TEST_ADD

#ifdef TEST_SCALAR
int main() {
	
	struct matrix *A = NULL;
	struct matrix *RES = NULL;
	struct matrix *REALRES = NULL;
	printf("Here \n");
	INIT_MATRIX(A, 3, 3);
	printf("After A \n");
	INIT_MATRIX(RES, 3, 3);
	printf("After RES \n");

	A->mat[0][0] = 6;
	A->mat[0][1] = 3;
	A->mat[0][2] = 0;
	
	A->mat[1][0] = 2;
	A->mat[1][1] = 5;
	A->mat[1][2] = 1;
	
	A->mat[2][0] = 9;
	A->mat[2][1] = 8;
	A->mat[2][2] = 6;
	
	printf("Printing A \n");
	printMatrix(A);
	printf("\n");
	
	int returnVal = scaleMatrix(A, 2, RES);
	
	if(returnVal != 0) {
		//print something
		printf("Failed with a code %d\n", returnVal);
	}
	else {
		//print the resulting matrix
		printf("Printing addition result \n");
		printMatrix(RES);
		printf("\n");
		INIT_MATRIX(REALRES, 3, 3);
		REALRES->mat[0][0] = 12;
		REALRES->mat[0][1] = 6;
		REALRES->mat[0][2] = 0;
		
		REALRES->mat[1][0] = 4;
		REALRES->mat[1][1] = 10;
		REALRES->mat[1][2] = 1;
		
		REALRES->mat[2][0] = 18;
		REALRES->mat[2][1] = 16;
		REALRES->mat[2][2] = 12;
		printf("Printing right addition result \n");
		printMatrix(REALRES);
		printf("\n");
		printf("Comparing results... and it %s \n", compareMatrieces(REALRES,RES)==0?"matches":"does not match");
	}
	assert(returnVal == 0);
	FREE_MATRIX(A);
	FREE_MATRIX(RES);
	FREE_MATRIX(REALRES);
	return 0;
}
#endif //TEST_SCALAR