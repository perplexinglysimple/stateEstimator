/**
 * @file EKF.h
 * @brief Extended Kalman Filter
 * 
 * This file contains the Extended Kalman Filter (EKF) algorithm.
 * 
 * This implents the Init, Predict, and Update functions for the EKF.
 * 
 * Note1: If the inverse operation produces NaN values then consider the following:
 *  1. Check for measurement outliers: Verify if any of the measurements used in the algorithm are outliers or 
 *     contain erroneous data. Outliers can significantly affect the covariance matrix and result in numerical 
 *     instability. Consider applying outlier rejection techniques or sensor fusion methods to mitigate the impact
 *     of outliers.
 *  2. Evaluate the system and measurement models: Review the equations and models used in the EKF algorithm. Ensure 
 *     that the system and measurement models are correctly formulated and represent the dynamics of the system accurately.
 *     A poorly defined or inaccurate model can lead to numerical instability and singularities.
 *  3. Implement numerical stability enhancements: Apply techniques to improve numerical stability, such as regularizing 
 *     the covariance matrix or employing numerical conditioning methods like the pseudo-inverse. These methods can help
 *     mitigate singularities and improve the stability of the algorithm.
 *  4. Adjust the process noise covariance: The process noise covariance matrix (Q) in the EKF algorithm represents the 
 *     uncertainty or variability in the system dynamics. If the algorithm is consistently encountering singularities, you
 *     may need to adjust the values in the process noise covariance matrix to ensure a well-conditioned system.
 *  5. Consider alternative algorithms: If the above steps do not resolve the issue, you may need to explore alternative 
 *     estimation algorithms. The EKF is just one approach, and there are other variants like the Unscented Kalman Filter (UKF) 
 *     or Particle Filters that can be more robust in handling non-linear or singular scenarios. I may implement these in the future.
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

/**
 * @brief Initialize EKF state from configuration.
 * @note If useFiniteDifferenceJacobian is false, jacobianF and jacobianH must be provided.
 * @note Set options.numberOfMeasurements to the measurement dimension (m). If 0, it defaults to n.
 */
EKFReturnCodes EKFInit(EKFState *ekf, EKFConfigOptions *options);

/**
 * @brief Predict step.
 * @note State propagation uses f(x). updateAMatrix/A are used for linear models
 *       and to support Jacobian computation.
 */
EKFReturnCodes EKFPredict(EKFState *ekf, double time, void* userData);

/**
 * @brief Update step with measurement.
 * @note Uses Joseph form covariance update for numerical stability.
 * @note Innovation matrix S inversion uses diagonal jitter if near-singular.
 */
EKFReturnCodes EKFUpdate(EKFState *ekf, EKFMeasurement *measurement);

/**
 * @brief Cleanup EKF state allocations (if mallocFlag is true).
 */
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
#define STATIC_ALLOC_EKF_DIRECTIVE(ekfptr, numstates, nummeas) \
    STATIC_MATRIX_DIRECTIVE(ekfptr->x, numstates, 1, x); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->P, numstates, numstates, P); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->Q, numstates, numstates, Q); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->R, nummeas, nummeas, R); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->A, numstates, numstates, A); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_P, numstates, numstates, _P); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_K, numstates, nummeas, _K); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_z, nummeas, 1, _z); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_F, numstates, numstates, _F); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_H, nummeas, numstates, _H); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_F_TRANSPOSE, numstates, numstates, _F_TRANSPOSE); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_H_TRANSPOSE, numstates, nummeas, _H_TRANSPOSE); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_x_predicted, numstates, 1, _x_predicted); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_S, nummeas, nummeas, _S); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_I, numstates, numstates, _I); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP1, numstates, numstates, _TEMP1); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP2, numstates, numstates, _TEMP2); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP3, numstates, 1, _TEMP3); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP4, nummeas, numstates, _TEMP4); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP5, numstates, nummeas, _TEMP5); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP6, nummeas, nummeas, _TEMP6); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP7, nummeas, 1, _TEMP7); \
    STATIC_MATRIX_DIRECTIVE(ekfptr->_TEMP8, numstates, 1, _TEMP8)

// This is used instead of the STATIC_ALLOC_EKF_DIRECTIVE macro and only allocates the non-_ prefixed matrices.
#define PRE_INIT_ALLOC(ekfptr, numstates, nummeas, malloc) \
    if (malloc) \
    { \
        INIT_MATRIX((ekfptr)->x, numstates, 1); \
        INIT_MATRIX((ekfptr)->P, numstates, numstates); \
        INIT_MATRIX((ekfptr)->Q, numstates, numstates); \
        INIT_MATRIX((ekfptr)->R, nummeas, nummeas); \
        INIT_MATRIX((ekfptr)->A, numstates, numstates); \
    } \
    else \
    { \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->x, numstates, 1, x); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->P, numstates, numstates, P); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->Q, numstates, numstates, Q); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->R, nummeas, nummeas, R); \
        STATIC_MATRIX_DIRECTIVE((ekfptr)->A, numstates, numstates, A); \
    }
    
#endif
