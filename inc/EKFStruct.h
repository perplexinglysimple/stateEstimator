#ifndef EKF_STRUCT_H
#define EKF_STRUCT_H

#include <stdbool.h>
#include <stdint.h>

#include "matrixMath.h"

// TODO fix matrix struct to not need this extra typedef
typedef struct matrix EKFMatrix;

typedef struct EKFState_ EKFState;

// State transition function prototype
typedef void (*EKFStateAFunction)(EKFMatrix* A, EKFMatrix* x, EKFState* ekf, double time, void* userData);
// System dynamics model function prototype
typedef void (*EKFStateTransitionFunction)(EKFMatrix* x, EKFMatrix* x_predicted, EKFState* ekf, void* userData);
// Measurement function prototype
typedef void (*EKFMeasurementFunction)(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
// Jacobian function prototype (fills J based on state x)
typedef void (*EKFJacobianFunction)(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);

struct EKFState_
{
    /** State vector (n x 1). */
    EKFMatrix* x;
    /** Predicted state vector (n x 1) produced by f(x). */
    EKFMatrix* _x_predicted;
    /** State covariance matrix (n x n). */
    EKFMatrix* P;
    /** Process noise covariance matrix (n x n). */
    EKFMatrix* Q;
    /** Measurement noise covariance matrix (m x m). */
    EKFMatrix* R;
    /** State transition matrix (n x n) for linear/locally linear models. */
    EKFMatrix* A;
    /** Temporary covariance storage (n x n). */
    EKFMatrix* _P;
    /** Kalman gain (n x m). */
    EKFMatrix* _K;
    /** Measurement residual vector storage (m x 1). */
    EKFMatrix* _z;
    /** Jacobian of the state transition function f (n x n). */
    EKFMatrix* _F;
    /** Jacobian of the measurement function h (m x n). */
    EKFMatrix* _H;
    /** Transpose of F (n x n). */
    EKFMatrix* _F_TRANSPOSE;
    /** Transpose of H (n x m). */
    EKFMatrix* _H_TRANSPOSE;
    /** Innovation covariance matrix S (m x m). */
    EKFMatrix* _S;
    /** Identity matrix (n x n). */
    EKFMatrix* _I;
    /** Temporary workspace matrix (n x n). */
    EKFMatrix* _TEMP1;
    /** Temporary workspace matrix (n x n). */
    EKFMatrix* _TEMP2;
    /** Temporary workspace vector (n x 1). */
    EKFMatrix* _TEMP3;
    /** Temporary workspace matrix (m x n). */
    EKFMatrix* _TEMP4;
    /** Temporary workspace matrix (n x m). */
    EKFMatrix* _TEMP5;
    /** Temporary workspace matrix (m x m). */
    EKFMatrix* _TEMP6;
    /** Temporary workspace vector (m x 1). */
    EKFMatrix* _TEMP7;
    /** Temporary workspace vector (n x 1). */
    EKFMatrix* _TEMP8;
    /** Callback to update A (used for linear models or as helper for Jacobians). */
    EKFStateAFunction updateAMatrix;
    /** State transition function f(x). */
    EKFStateTransitionFunction f;
    /** Measurement function h(x). */
    EKFMeasurementFunction h;
    /** Jacobian of f(x) with respect to x (optional if finite-diff enabled). */
    EKFJacobianFunction jacobianF;
    /** Jacobian of h(x) with respect to x (optional if finite-diff enabled). */
    EKFJacobianFunction jacobianH;
    /** Number of state variables (n). */
    int numberOfStates;
    /** Number of measurement variables (m). */
    int numberOfMeasurements;
    /** Use finite-difference Jacobians when true; otherwise jacobianF/jacobianH must be provided. */
    bool useFiniteDifferenceJacobian;
    /** True if EKF matrices were allocated with malloc and must be freed in cleanup. */
    bool mallocFlag;
};

typedef struct EKFConfigOptions_
{
    /** Initial state vector x0 (n x 1). */
    EKFMatrix* x0;
    /** Initial covariance matrix P0 (n x n). */
    EKFMatrix* P0;
    /** Process noise covariance matrix Q (n x n). */
    EKFMatrix* Q;
    /** Measurement noise covariance matrix R (m x m). */
    EKFMatrix* R;
    /** Initial/constant state transition matrix A (n x n). */
    EKFMatrix* A;
    /** Number of state variables (n). */
    int n;
    /** Callback to update A for linear models or as a helper for Jacobians. */
    EKFStateAFunction updateAMatrix;
    /** State transition function f(x). */
    EKFStateTransitionFunction f;
    /** Measurement function h(x). */
    EKFMeasurementFunction h;
    /** Jacobian of f(x) with respect to x (required if finite-diff disabled). */
    EKFJacobianFunction jacobianF;
    /** Jacobian of h(x) with respect to x (required if finite-diff disabled). */
    EKFJacobianFunction jacobianH;
    /** Optional alias of n (kept for backwards compatibility). */
    int numberOfStates;
    /** Measurement dimension m (if 0, defaults to n). */
    int numberOfMeasurements;
    /** Use finite-difference Jacobians when true; otherwise jacobianF/jacobianH must be provided. */
    bool useFiniteDifferenceJacobian;
    /** True to malloc internal matrices; false to use static matrices. */
    bool mallocFlag;
} EKFConfigOptions;

typedef struct EKFMeasurement_
{
    /** Measurement vector z (m x 1). */
    EKFMatrix* z;
} EKFMeasurement;

#endif // EKF_STRUCT_H
