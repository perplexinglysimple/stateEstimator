/**
 * @file matrixMath.h
 * @brief Matrix math types, macros, and operations used by the EKF.
 */
#ifndef MATRIX_MATH
#define MATRIX_MATH

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>

#include "utils.h"

/** @brief Numeric precision used throughout the EKF math layer. */
// CHANGE THIS TO FLOAT OR DOUBLE DEPENDING ON YOUR NEEDS
typedef double ekfType;

// Enable extra matrix math checks in debug builds (disabled in release).
// Define MATRIX_MATH_DEBUG_CHECKS=1 to force-enable in any build.
#ifndef MATRIX_MATH_DEBUG_CHECKS
#if defined(NDEBUG)
#define MATRIX_MATH_DEBUG_CHECKS 0
#else
#define MATRIX_MATH_DEBUG_CHECKS 1
#endif
#endif

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

/** @brief Return codes for matrix operations. */
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

/**
 * @brief Dense matrix container.
 *
 * Matrices can be allocated in "jagged" mode (array of row pointers) or
 * in a contiguous static buffer via STATIC_MATRIX_DIRECTIVE.
 */
struct matrix {
	matrixType **mat;
	matrixType *_mat;
	int initilized;
	int row;
	int col;
	bool jaggedAlloc;
};

/** @brief Allocate and zero-initialize a jagged matrix. */
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

#define MM_CONCAT(a, b) MM_CONCAT_INNER(a, b)
#define MM_CONCAT_INNER(a, b) a ## b

#define MM_UNIQUE_NAME(base) MM_CONCAT(base, __COUNTER__)
#define UNIQUE_NAME_PER_MACRO(base) MM_CONCAT(base, __LINE__)

/** @brief Access element (i, j) from a static matrix. */
#define ACCESS_STATIC_MATRIX(m, i, j) \
	(m)._mat[i * (m).col + j]

/** @brief Access element (i, j) from a matrix, regardless of allocation mode. */
#define ACCESS_MATRIX(m, i, j) \
	(matrixType) ((m).jaggedAlloc ? (m).mat[i][j] : ACCESS_STATIC_MATRIX(m, i, j))

/** @brief Set element (i, j) in a matrix, regardless of allocation mode. */
#define SET_MATRIX(m, i, j, val) \
	do { \
		if ((m).jaggedAlloc) { \
			(m).mat[i][j] = (val); \
		} else { \
			ACCESS_STATIC_MATRIX((m), i, j) = (val); \
		} \
	} while (0)

#define ZEROIZE_MATRIX(m, i, j) \
	for(i = 0; i < (m).row; ++i) { \
		for(j = 0; j < (m).col; ++j) { \
			SET_MATRIX((m), i, j, 0); \
		} \
	}

#define _ZEROIZE_MATRIX(m, i, j) \
	for(i = 0; i < m.row; ++i) { \
		for(j = 0; j < m.col; ++j) { \
			ACCESS_STATIC_MATRIX(m, i, j) = 0; \
		} \
	}

/** @brief Copy a 2D C array into a matrix. */
#define COPY_2DARRAY_TO_MATRIX(array, matrixptr) \
	do { \
		int _i, _j; \
		for(_i = 0; _i < matrixptr->row; ++_i) { \
			for(_j = 0; _j < matrixptr->col; ++_j) { \
				SET_MATRIX((*matrixptr), _i, _j, array[_i][_j]); \
			} \
		} \
	} while(0)

/**
 * @brief Create a static matrix with storage backing.
 *
 * This macro declares a static 2D array and binds a struct matrix to it.
 * It is not thread-safe and should not be used in concurrent contexts.
 */
#define STATIC_MATRIX_DIRECTIVE(ptr, _row, _col, name) \
	static matrixType MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat))[(_row)][(_col)]; \
	static struct matrix MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m)) = { .mat = NULL, .initilized = 1, .row = (_row), .col = (_col), .jaggedAlloc = false }; \
	MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat = (matrixType **)MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m_mat)); \
	MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m))._mat = (matrixType *)MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m)).mat; \
	int MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(i)), MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(j)); \
	_ZEROIZE_MATRIX(MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m)), MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(i)), MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(j))); \
	ptr = &MM_CONCAT(name, UNIQUE_NAME_PER_MACRO(m));
	

/** @brief Free a jagged matrix allocated with INIT_MATRIX. */
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

#if MATRIX_MATH_DEBUG_CHECKS
#define NULL_CHECK_MATRIX(ptr) \
	do { \
		if(ptr == NULL) { \
			LOG_ERROR("NULL_CHECK_MATRIX: ptr is NULL"); \
			return MATRIX_NULL_POINTER; \
		} \
	} while(0)
#else
#define NULL_CHECK_MATRIX(ptr)
#endif

#if MATRIX_MATH_DEBUG_CHECKS
#define NULL_CHECK_MATRIX_RES(ptr) \
	do { \
		if(ptr == NULL) { \
			LOG_ERROR("NULL_CHECK_MATRIX_RES: ptr is NULL"); \
			return MATRIX_NULL_RES_POINTER; \
		} \
	} while(0)
#else
#define NULL_CHECK_MATRIX_RES(ptr)
#endif

#if MATRIX_MATH_DEBUG_CHECKS
#define NON_INIT_CHECK_MATRIX(ptr) \
	do { \
		if(ptr->initilized == 0) { \
			LOG_ERROR("NON_INIT_CHECK_MATRIX: ptr is not initilized"); \
			return MATRIX_NOT_INITIALIZED; \
		} \
	} while(0)
#else
#define NON_INIT_CHECK_MATRIX(ptr)
#endif

// Guard against in-place usage unless explicitly allowed.
#if MATRIX_MATH_DEBUG_CHECKS
#define NO_ALIAS_CHECK_MATRIX2(a, res) \
	do { \
		if ((res) == (a)) { \
			LOG_ERROR("NO_ALIAS_CHECK_MATRIX2: res must not alias input"); \
			return MATRIX_ERROR; \
		} \
	} while(0)

#define NO_ALIAS_CHECK_MATRIX3(a, b, res) \
	do { \
		if ((res) == (a) || (res) == (b)) { \
			LOG_ERROR("NO_ALIAS_CHECK_MATRIX3: res must not alias inputs"); \
			return MATRIX_ERROR; \
		} \
	} while(0)
#else
#define NO_ALIAS_CHECK_MATRIX2(a, res)
#define NO_ALIAS_CHECK_MATRIX3(a, b, res)
#endif

#if MATRIX_MATH_DEBUG_CHECKS
#define DIMENSION_CHECK_MULT_MATRIX(a, b, res) \
	do { \
		if(a->col != b->row || a->row != res->row || b->col != res->col) { \
			LOG_ERROR("DIMENSION_CHECK_MULT_MATRIX: a->col != b->row || a->row != res->row || b->col != res->col"); \
			LOG_INFO("a->col = %d, b->row = %d, a->row = %d, res->row = %d, b->col = %d, res->col = %d", a->col, b->row, a->row, res->row, b->col, res->col); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)
#else
#define DIMENSION_CHECK_MULT_MATRIX(a, b, res)
#endif

#if MATRIX_MATH_DEBUG_CHECKS
#define DIMENSION_CHECK_SCALER_MATRIX(a, res) \
	do { \
		if(a->col != res->col || a->row != res->row) { \
			LOG_ERROR("DIMENSION_CHECK_SCALER_MATRIX: a->col != res->col || a->row != res->row"); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)
#else
#define DIMENSION_CHECK_SCALER_MATRIX(a, res)
#endif
		
#if MATRIX_MATH_DEBUG_CHECKS
#define DIMENSION_CHECK_ADD_MATRIX(a, b, res) \
	do { \
		if(a->col != b->col || a->row != b->row || a->col != res->col || a->row != res->row) { \
			LOG_ERROR("DIMENSION_CHECK_ADD_MATRIX: a->col != b->col || a->row != b->row || a->col != res->col || a->row != res->row"); \
			return MATRIX_DIMENSION_MISMATCH; \
		} \
	} while(0)
#else
#define DIMENSION_CHECK_ADD_MATRIX(a, b, res)
#endif

#if MATRIX_MATH_DEBUG_CHECKS
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
/**
 * @brief Multiply two matrices (res = a * b).
 * @note Not safe for in-place usage. res must not alias a or b.
 */
matrixReturnCodes multMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

/**
 * @brief Scale a matrix (res = a * scaler).
 * @note Safe for in-place usage when res == a.
 */
matrixReturnCodes scaleMatrix(struct matrix *a, struct matrix *res, matrixType scaler);

/**
 * @brief Add two matrices (res = a + b).
 * @note Not safe for in-place usage. res must not alias a or b.
 */
matrixReturnCodes addMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

/**
 * @brief Subtract two matrices (res = a - b).
 * @note Not safe for in-place usage. res must not alias a or b.
 */
matrixReturnCodes subMatrix(struct matrix *a, struct matrix *b, struct matrix *res);

/**
 * @brief Invert a matrix (res = a^-1).
 * @note Not safe for in-place usage. res must not alias a.
 */
matrixReturnCodes inverseMatrix(struct matrix *a, struct matrix *res);

/**
 * @brief Invert a matrix with diagonal jitter retry (res = a^-1).
 * @note Not safe for in-place usage. res must not alias a.
 * @note This function mutates a by adding jitter to the diagonal on failed attempts.
 */
matrixReturnCodes inverseMatrixWithJitter(struct matrix *a, struct matrix *res, matrixType jitter, int maxAttempts, matrixType jitterScale);

/**
 * @brief Set matrix to identity.
 * @note Only valid for square matrices.
 */
matrixReturnCodes setIdentityMatrix(struct matrix *a);

/**
 * @brief Transpose a matrix (res = a^T).
 * @note Not safe for in-place usage. res must not alias a.
 */
matrixReturnCodes transposeMatrix(struct matrix *a, struct matrix *b);

/**
 * @brief Compute (I - A) into res.
 * @note Not safe for in-place usage. res must not alias a.
 */
matrixReturnCodes identityMatrixMinusA(struct matrix *a, struct matrix *res);

/**
 * @brief Compare two matrices for exact equality.
 * @note Safe when a and b alias.
 */
matrixReturnCodes compareMatrieces(struct matrix *a, struct matrix *b);

/**
 * @brief Print a matrix to stdout.
 */
void printMatrix(struct matrix *a);

/**
 * @brief Copy matrix a into res.
 * @note Not safe for in-place usage. res must not alias a.
 */
matrixReturnCodes copyMatrix(struct matrix *a, struct matrix *res);

/**
 * @brief Check matrix for NaN values.
 */
matrixReturnCodes nanCheckMatrix(struct matrix *a);

#endif
