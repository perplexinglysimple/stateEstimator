#ifndef MATRIX_MATH
#define MATRIX_MATH

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "utils.h"

#ifndef MATRIXTYPE
#define MATRIXTYPE
typedef int matrixType;
//This is the type so the print statements wont bug out
//This can be int or double
#define _INT
#endif

typedef enum matrixReturnCodes_ {
	MATRIX_SUCCESS = 0,
	MATRIX_ERROR = -1,
	MATRIX_NULL_POINTER = -2,
	MATRIX_NULL_RES_POINTER = -3,
	MATRIX_DIMENSION_MISMATCH = -4,
	MATRIX_NOT_INITIALIZED = -5,
	MATRIX_COMPARE_FAILURE = -6,
} matrixReturnCodes;

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

#define __CONCAT(a, b) __CONCAT_INNER(a, b)
#define __CONCAT_INNER(a, b) a ## b

#define UNIQUE_NAME(base) __CONCAT(base, __COUNTER__)
#define UNIQUE_NAME_PER_MACRO(base) __CONCAT(base, __LINE__)

#define STATIC_MATRIX_DIRECTIVE(ptr, _row, _col, name) \
	static matrixType __CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat))[_row][_col]; \
	static struct matrix __CONCAT(name, UNIQUE_NAME_PER_MACRO(m)) = { .mat = NULL, .initilized = 1, .row = _row, .col = _col }; \
	__CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat = (matrixType **)__CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat)); \
	int __CONCAT(name, UNIQUE_NAME_PER_MACRO(i)), __CONCAT(name, UNIQUE_NAME_PER_MACRO(j)); \
	for(__CONCAT(name, UNIQUE_NAME_PER_MACRO(i)) = 0; __CONCAT(name, UNIQUE_NAME_PER_MACRO(i)) < _row; ++__CONCAT(name, UNIQUE_NAME_PER_MACRO(i))) { \
		for(__CONCAT(name, UNIQUE_NAME_PER_MACRO(j)) = 0; __CONCAT(name, UNIQUE_NAME_PER_MACRO(j)) < _col; ++__CONCAT(name, UNIQUE_NAME_PER_MACRO(j))) { \
			__CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat[__CONCAT(name, UNIQUE_NAME_PER_MACRO(i))][__CONCAT(name, UNIQUE_NAME_PER_MACRO(j))] = 0; \
		} \
	} \
	ptr = &__CONCAT(name, UNIQUE_NAME_PER_MACRO(m));
	
		
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

#define NULL_CHECK_MATRIX(ptr) \
	do { \
		if(ptr == NULL) { \
			LOG_ERROR("NULL_CHECK_MATRIX: ptr is NULL"); \
			return MATRIX_NULL_POINTER; \
		} \
	} while(0)

#define NULL_CHECK_MATRIX_RES(ptr) \
	do { \
		if(ptr == NULL) { \
			LOG_ERROR("NULL_CHECK_MATRIX_RES: ptr is NULL"); \
			return MATRIX_NULL_RES_POINTER; \
		} \
	} while(0)

#define NON_INIT_CHECK_MATRIX(ptr) \
	do { \
		if(ptr->initilized == 0) { \
			LOG_ERROR("NON_INIT_CHECK_MATRIX: ptr is not initilized"); \
			return MATRIX_NOT_INITIALIZED; \
		} \
	} while(0)

#define DIMENSION_CHECK_MULT_MATRIX(a, b, res) \
	do { \
		if(a->col != b->row || a->row != res->row || b->col != res->col) { \
			LOG_ERROR("DIMENSION_CHECK_MULT_MATRIX: a->col != b->row || a->row != res->row || b->col != res->col"); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)

#define DIMENSION_CHECK_SCALER_MATRIX(a, res) \
	do { \
		if(a->col != res->col || a->row != res->row) { \
			LOG_ERROR("DIMENSION_CHECK_SCALER_MATRIX: a->col != res->col || a->row != res->row"); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)
		
#define DIMENSION_CHECK_ADD_MATRIX(a, b, res) \
	do { \
		if(a->col != b->col || a->row != b->row || a->col != res->col || a->row != res->row) { \
			LOG_ERROR("DIMENSION_CHECK_ADD_MATRIX: a->col != b->col || a->row != b->row || a->col != res->col || a->row != res->row"); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)
	
	
// ------------------------- Public Functions ------------------------- //
matrixReturnCodes multMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes scaleMatrix(struct matrix *a, struct matrix *res, matrixType scaler);

matrixReturnCodes addMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes subMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes inverseMatrix(struct matrix *a, struct matrix *res);

matrixReturnCodes tranposeMatrix(struct matrix *a, struct matrix *b);

matrixReturnCodes idenityMatrixMinusA(struct matrix *a, struct matrix *res);

matrixReturnCodes compareMatrieces(struct matrix *a, struct matrix *b);

void printMatrix(struct matrix *a);

matrixReturnCodes copyMatrix(struct matrix *a, struct matrix *res);

#endif