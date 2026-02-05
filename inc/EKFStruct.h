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

typedef struct EKFState_ {
    EKFMatrix* x; // State vector
    EKFMatrix* _x_predicted; // State vector predicted by the state transition function
    EKFMatrix* P; // Covariance matrix
    EKFMatrix* Q; // Process noise covariance matrix
    EKFMatrix* R; // Measurement noise covariance matrix
    EKFMatrix* A; // State transition matrix
    EKFMatrix* _P; // Covariance matrix temporary storage
    EKFMatrix* _K; // Kalman gain
    EKFMatrix* _z; // State transition matrix
    EKFMatrix* _F; // Jacobian of the state transition function
    EKFMatrix* _H; // Jacobian of the measurement function
    EKFMatrix* _F_TRANSPOSE; // Jacobian of the state transition function transposed
    EKFMatrix* _H_TRANSPOSE; // Jacobian of the measurement function transposed
    EKFMatrix* _S; // Innovation covariance matrix
    EKFMatrix* _I; // Identity matrix
    EKFMatrix* _TEMP1; // Temporary matrix
    EKFMatrix* _TEMP2; // Temporary matrix
    EKFMatrix* _TEMP3; // Temporary vector
    EKFMatrix* _TEMP4; // Temporary matrix (m x n)
    EKFMatrix* _TEMP5; // Temporary matrix (n x m)
    EKFMatrix* _TEMP6; // Temporary matrix (m x m)
    EKFMatrix* _TEMP7; // Temporary vector (m x 1)
    EKFMatrix* _TEMP8; // Temporary vector (n x 1)
    EKFStateAFunction updateAMatrix; // State transition matrix function
    EKFStateTransitionFunction f; // System dynamics model function
    EKFMeasurementFunction h; // Measurement function
    EKFJacobianFunction jacobianF; // Jacobian of the state transition function
    EKFJacobianFunction jacobianH; // Jacobian of the measurement function
    int numberOfStates; // Number of state variables
    int numberOfMeasurements; // Number of measurement variables
    bool useFiniteDifferenceJacobian; // Use finite difference to calculate the Jacobian. If false, jacobianF/jacobianH must be provided.
    bool mallocFlag; // Flag to indicate if the matrices were malloced or not
} EKFState;

typedef struct EKFConfigOptions_ {
    EKFMatrix* x0; // Initial state vector
    EKFMatrix* P0; // Initial covariance matrix
    EKFMatrix* Q; // Process noise covariance matrix
    EKFMatrix* R; // Measurement noise covariance matrix
    EKFMatrix* A; // State transition matrix
    int n; // Number of state variables
    EKFStateAFunction updateAMatrix; // State transition matrix function
    EKFStateTransitionFunction f; // System dynamics model function
    EKFMeasurementFunction h; // Measurement function
    EKFJacobianFunction jacobianF; // Jacobian of the state transition function
    EKFJacobianFunction jacobianH; // Jacobian of the measurement function
    int numberOfStates; // Number of state variables
    int numberOfMeasurements; // Number of measurement variables
    bool useFiniteDifferenceJacobian; // Use finite difference to calculate the Jacobian. If false, jacobianF/jacobianH must be provided.
    bool mallocFlag; // Flag to indicate if the matrices were malloced or not
} EKFConfigOptions;

typedef struct EKFMeasurement_ {
    EKFMatrix* z; // Measurement vector
} EKFMeasurement;

#endif // EKF_STRUCT_H
