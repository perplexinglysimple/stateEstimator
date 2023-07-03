#ifndef MATRIX_MATH
#define MATRIX_MATH

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include "utils.h"

// CHANGE THIS TO FLOAT OR DOUBLE DEPENDING ON YOUR NEEDS
typedef double ekfType;

// If defined, the matrix math functions will check for NaNs and infinities
// This is useful for debugging, but will slow down the code
#define NAN_CHECK

#ifndef MATRIXTYPE
#define MATRIXTYPE
typedef ekfType matrixType;
#define _DOUBLE
#endif

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
	MATRIX_NAN_FAILURE = -7,
	MATRIX_INF_FAILURE = -8,
} matrixReturnCodes;

struct matrix {
	matrixType **mat;
	matrixType *_mat;
	int initilized;
	int row;
	int col;
	bool jaggedAlloc;
};

#define INIT_MATRIX(ptr, _row, _col) \
	do { \
	    FREE_MATRIX(ptr); \
	    ptr = malloc(1 * sizeof(struct matrix)); \
		ptr->col = _col; \
		ptr->row = _row; \
		ptr->initilized = 1; \
		ptr->mat = NULL; \
		ptr->jaggedAlloc = true; \
		int _i, _j; \
		ptr->mat = malloc(_row * sizeof(matrixType *)); \
		for(_i = 0; _i < _row; ++_i) { \
			ptr->mat[_i] = malloc(_col * sizeof(matrixType)); \
			for(_j = 0; _j < _col; ++_j) { \
				ptr->mat[_i][_j] = 0; \
			} \
		} \
	} while(0)

#define __CONCAT(a, b) __CONCAT_INNER(a, b)
#define __CONCAT_INNER(a, b) a ## b

#define UNIQUE_NAME(base) __CONCAT(base, __COUNTER__)
#define UNIQUE_NAME_PER_MACRO(base) __CONCAT(base, __LINE__)

#define ACCESS_STATIC_MATRIX(m, i, j) \
	(m)._mat[i * (m).row + j]

#define ACCESS_MATRIX(m, i, j) \
	(matrixType) ((m).jaggedAlloc ? (m).mat[i][j] : ACCESS_STATIC_MATRIX(m, i, j))

#define ZEROIZE_MATRIX(m, i, j, _row, _col) \
	for(i = 0; i < m.row; ++i) { \
		for(j = 0; j < m.col; ++j) { \
			ACCESS_MATRIX(m, i, j) = 0; \
		} \
	}

#define _ZEROIZE_MATRIX(m, i, j, _row, _col) \
	for(i = 0; i < m.row; ++i) { \
		for(j = 0; j < m.col; ++j) { \
			ACCESS_STATIC_MATRIX(m, i, j) = 0; \
		} \
	}

#define STATIC_MATRIX_DIRECTIVE(ptr, _row, _col, name) \
	static matrixType __CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat))[(_row)][(_col)]; \
	static struct matrix __CONCAT(name, UNIQUE_NAME_PER_MACRO(m)) = { .mat = NULL, .initilized = 1, .row = (_row), .col = (_col), .jaggedAlloc = false }; \
	__CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat = (matrixType **)__CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat)); \
	__CONCAT(name, UNIQUE_NAME_PER_MACRO(m))._mat = (matrixType *)__CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat; \
	int __CONCAT(name, UNIQUE_NAME_PER_MACRO(i)), __CONCAT(name, UNIQUE_NAME_PER_MACRO(j)); \
	_ZEROIZE_MATRIX(__CONCAT(name, UNIQUE_NAME_PER_MACRO(m)), __CONCAT(name, UNIQUE_NAME_PER_MACRO(i)), __CONCAT(name, UNIQUE_NAME_PER_MACRO(j)), (_row), (_col)); \
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
			LOG_INFO("a->col = %d, b->row = %d, a->row = %d, res->row = %d, b->col = %d, res->col = %d", a->col, b->row, a->row, res->row, b->col, res->col); \
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

#if defined(NAN_CHECK)
#define NAN_CHECK_MATRIX(a) \
	do { \
		if(nanCheckMatrix(a) != MATRIX_SUCCESS) { \
			LOG_ERROR("NAN_CHECK_MATRIX: nanCheckMatrix(a) != MATRIX_SUCCESS"); \
			return MATRIX_NAN_FAILURE; \
		} \
	} while(0)
#else
#define NAN_CHECK_MATRIX(a)
#endif
	
	
// ------------------------- Public Functions ------------------------- //
matrixReturnCodes multMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes scaleMatrix(struct matrix *a, struct matrix *res, matrixType scaler);

matrixReturnCodes addMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes subMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

matrixReturnCodes inverseMatrix(struct matrix *a, struct matrix *res);

matrixReturnCodes transposeMatrix(struct matrix *a, struct matrix *b);

matrixReturnCodes identityMatrixMinusA(struct matrix *a, struct matrix *res);

matrixReturnCodes compareMatrieces(struct matrix *a, struct matrix *b);

void printMatrix(struct matrix *a);

matrixReturnCodes copyMatrix(struct matrix *a, struct matrix *res);

matrixReturnCodes nanCheckMatrix(struct matrix *a);

#endif