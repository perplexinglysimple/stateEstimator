#ifndef EKF_STRUCT_H
#define EKF_STRUCT_H

// CHANGE THIS TO FLOAT OR DOUBLE DEPENDING ON YOUR NEEDS
typedef double ekfType;

#define MATRIXTYPE
typedef ekfType matrixType;
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
} EKFState;

typedef struct EKFConfigOptions_ {
    EKFMatrix* x0; // Initial state vector
    EKFMatrix* P0; // Initial covariance matrix
    EKFMatrix* Q; // Process noise covariance matrix
    EKFMatrix* R; // Measurement noise covariance matrix
    EKFMatrix* A; // State transition matrix
    int n; // Number of state variables
} EKFConfigOptions;

typedef struct EKFMeasurement_ {
    EKFMatrix* z; // Measurement vector
    EKFMatrix* H; // Measurement matrix
} EKFMeasurement;

#endif // EKF_STRUCT_H