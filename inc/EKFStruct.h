#ifndef EKF_STRUCT_H
#define EKF_STRUCT_H

#include <stdbool.h>

#include "matrixMath.h"

// TODO fix matrix struct to not need this extra typedef
typedef struct matrix EKFMatrix;

// System dynamics model function prototype
typedef void (*EKFStateTransitionFunction)(EKFMatrix* x, EKFMatrix* x_predicted);
// Measurement function prototype
typedef void (*EKFMeasurementFunction)(EKFMatrix* x, EKFMatrix* z);

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
    EKFStateTransitionFunction f; // System dynamics model function
    EKFMeasurementFunction h; // Measurement function
    int numberOfStates; // Number of state variables
    bool useFiniteDifferenceJacobian; // Use finite difference to calculate the Jacobian. If false, the user must provide the Jacobian through the _F and _H matrices.
    bool mallocFlag; // Flag to indicate if the matrices were malloced or not
} EKFState;

typedef struct EKFConfigOptions_ {
    EKFMatrix* x0; // Initial state vector
    EKFMatrix* P0; // Initial covariance matrix
    EKFMatrix* Q; // Process noise covariance matrix
    EKFMatrix* R; // Measurement noise covariance matrix
    EKFMatrix* A; // State transition matrix
    int n; // Number of state variables
    EKFStateTransitionFunction f; // System dynamics model function
    EKFMeasurementFunction h; // Measurement function
    int numberOfStates; // Number of state variables
    bool useFiniteDifferenceJacobian; // Use finite difference to calculate the Jacobian. If false, the user must provide the Jacobian through the _F and _H matrices.
    bool mallocFlag; // Flag to indicate if the matrices were malloced or not
} EKFConfigOptions;

typedef struct EKFMeasurement_ {
    EKFMatrix* z; // Measurement vector
} EKFMeasurement;

#endif // EKF_STRUCT_H