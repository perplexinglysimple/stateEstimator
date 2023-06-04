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
#define EPSILON (0.0001)

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

#define CHECK_

#endif