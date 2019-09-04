#ifndef MATRIX_MATH
#define MATRIX_MATH

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#ifndef MATRIXTYPE
#define MATRIXTYPE
typedef int matrixType;
//This is the type so the print statements wont bug outp
//This can be int or double
#define _INT
#endif


struct matrix {
	matrixType **mat;
	int initilized;
	int row;
	int col;
};

#define INIT_MATRIX(ptr, _row, _col) \
	do { \
	    FREE_MATRIX(ptr); \
	    ptr = calloc(1, sizeof(struct matrix)); \
		ptr->col = _col; \
		ptr->row = _row; \
		ptr->initilized = 1; \
		ptr->mat = NULL; \
		int _i, _j; \
		ptr->mat = calloc(_row, sizeof(matrixType *)); \
		for(_i = 0; _i < _row; ++_i) { \
			ptr->mat[_i] = calloc(_col, sizeof(matrixType) ); \
			for(_j = 0; _j < _col; ++_j) { \
				ptr->mat[_i][_j] = 0; \
			} \
		} \
	} while(0)
		
#define FREE_MATRIX(ptr) \
	do { \
	    if(ptr != NULL) { \
		  int row = ptr->row; \
		  ptr->col = 0; \
		  ptr->row = 0; \
		  ptr->initilized = 0; \
		  int _i; \
		  for(_i = 0; _i < row; ++_i) { \
			free(ptr->mat[_i]); \
		  } \
		  free(ptr->mat); \
		  ptr->mat = NULL; \
		  free(ptr); \
		  ptr = NULL; \
		} \
	} while(0)
		
	
	
//Inputs are three matrixes. a and b matrix will be multiplied and the result will be output to res
int multMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

int scaleMatrix(struct matrix *a, struct matrix *res, matrixType scaler);

int addMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

void printMatrix(struct matrix *a);

int compareMatrieces(struct matrix *a, struct matrix *b);

int tranposeMatrix(struct matrix *a, struct matrix *b);

	

#endif