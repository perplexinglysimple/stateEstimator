/**
 * @file EKF.h
 * @brief Extended Kalman Filter
 * 
 * This file contains the Extended Kalman Filter (EKF) algorithm.
 * 
 * This implents the Init, Predict, and Update functions for the EKF.
 * 
*/

#ifndef EKF_H
#define EKF_H

#include "EKFStruct.h"
#include "matrixMath.h"
#include "utils.h"

// Epsilon is used for small perturbations for finite difference Jacobian calculations.
#define EPSILON ((matrixType) 0.0001)

typedef enum EKFReturnCodes_ {
    EKF_SUCCESS = 0,
    EKF_ERROR = -1,
    EKF_NULL_POINTER = -2,
} EKFReturnCodes;

EKFReturnCodes EKFInit(EKFState *ekf, EKFConfigOptions *options);
EKFReturnCodes EKFPredict(EKFState *ekf);
EKFReturnCodes EKFUpdate(EKFState *ekf, EKFMeasurement *measurement);
EKFReturnCodes EKFCleanup(EKFState *ekf);

#define NULL_CHECK_EKF(ekf) \
    do { \
        if (ekf == NULL) { \
            LOG_ERROR("Null pointer passed to EKF function."); \
            return EKF_NULL_POINTER; \
        } \
    } while (0)

#define MATRIX_MATH_RETURN_CHECK(ret) \
    do { \
        if (ret != MATRIX_SUCCESS) { \
            LOG_ERROR("Matrix math error in EKF function."); \
            return EKF_ERROR; \
        } \
    } while (0)

// This macro is meant as a helper for the EKFInit() function in tests
// but is not robust enough to be used outside of that context. Use at your own risk.
#define STATIC_ALLOC_EKF_DIRECTIVE(ekfptr, numstates) \
    STATIC_MATRIX_DIRECTIVE(ekfptr->x, numstates, 1, x); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->P, numstates, numstates, P); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->Q, numstates, numstates, Q); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->R, numstates, numstates, R); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->A, numstates, numstates, A); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_P, numstates, numstates, _P); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_K, numstates, numstates, _K); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_z, numstates, 1, _z); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_F, numstates, numstates, _F); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_H, numstates, numstates, _H); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_F_TRANSPOSE, numstates, numstates, _F_TRANSPOSE); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_H_TRANSPOSE, numstates, numstates, _H_TRANSPOSE)

// This is used instead of the STATIC_ALLOC_EKF_DIRECTIVE macro and only allocates the non-_ prefixed matrices.
#define PRE_INIT_ALLOC(ekfptr, numstates, malloc) \
    if (malloc) \
    { \
        INIT_MATRIX((ekfptr)->x, numstates, 1); \
        INIT_MATRIX((ekfptr)->P, numstates, numstates); \
        INIT_MATRIX((ekfptr)->Q, numstates, numstates); \
        INIT_MATRIX((ekfptr)->R, numstates, numstates); \
        INIT_MATRIX((ekfptr)->A, numstates, numstates); \
    } \
    else \
    { \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->x, numstates, 1, x); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->P, numstates, numstates, P); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->Q, numstates, numstates, Q); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->R, numstates, numstates, R); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->A, numstates, numstates, A); \
    }
    
#endif