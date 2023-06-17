#include "matrixMath.h"
#include "EKF.h"

// ------------------------- Private Function Prototypes ------------------------- //
/// @brief Calculate the Jacobian matrix, F, of the state transition function with respect to the state variables.
void calculateJacobian(EKFMatrix *x, EKFMatrix *_x_predicted, EKFMatrix *Jacobian, EKFStateTransitionFunction f, int numOfStates);

// ------------------------- Public Functions ------------------------- //
// See EKF.h for documentation
EKFReturnCodes EKFInit(EKFState *ekf, EKFConfigOptions *options)
{
  NULL_CHECK_EKF(ekf);
  NULL_CHECK_EKF(options);

  NULL_CHECK_MATRIX(options->x0);
  NULL_CHECK_MATRIX(options->P0);
  NULL_CHECK_MATRIX(options->Q);
  NULL_CHECK_MATRIX(options->R);
  NULL_CHECK_MATRIX(options->A);

  // Check to see if x, P, Q, R, and A are the correct size.
  if (options->x0->row != options->n || options->x0->col != 1)
  {
    LOG_ERROR("EKFInit() failed because the x0 matrix is not the correct size.");
    return EKF_ERROR;
  }
  if (options->P0->row != options->n || options->P0->col != options->n)
  {
    LOG_ERROR("EKFInit() failed because the P0 matrix is not the correct size.");
    return EKF_ERROR;
  }
  if (options->Q->row != options->n || options->Q->col != options->n)
  {
    LOG_ERROR("EKFInit() failed because the Q matrix is not the correct size.");
    return EKF_ERROR;
  }
  if (options->R->row != options->n || options->R->col != options->n)
  {
    LOG_ERROR("EKFInit() failed because the R matrix is not the correct size.");
    return EKF_ERROR;
  }
  if (options->A->row != options->n || options->A->col != options->n)
  {
    LOG_ERROR("EKFInit() failed because the A matrix is not the correct size.");
    return EKF_ERROR;
  }
  // Check to see if f and h are not null.
  NULL_CHECK_EKF(options->f);
  NULL_CHECK_EKF(options->h);

  // Initialize the state vector.
  MATRIX_MATH_RETURN_CHECK(copyMatrix(options->x0, ekf->x));
  // Initialize the state covariance matrix.
  MATRIX_MATH_RETURN_CHECK(copyMatrix(options->P0, ekf->P));
  // Initialize the process noise covariance matrix.
  MATRIX_MATH_RETURN_CHECK(copyMatrix(options->Q, ekf->Q));
  // Initialize the measurement noise covariance matrix.
  MATRIX_MATH_RETURN_CHECK(copyMatrix(options->R, ekf->R));
  // Initialize the state transition matrix.
  MATRIX_MATH_RETURN_CHECK(copyMatrix(options->A, ekf->A));
  // Initialize the temporary storage matrices.
  ekf->mallocFlag = options->mallocFlag;
  if (ekf->mallocFlag)
  {
    INIT_MATRIX(ekf->_P, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_K, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_z, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_F, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_H, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_F_TRANSPOSE, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_H_TRANSPOSE, ekf->P->row, ekf->P->col);
  }
  else
  {
    // Check that the temporary storage matrices are the correct size and not null
    NULL_CHECK_MATRIX(ekf->_P);
    NULL_CHECK_MATRIX(ekf->_K);
    NULL_CHECK_MATRIX(ekf->_z);
    NULL_CHECK_MATRIX(ekf->_F);
    NULL_CHECK_MATRIX(ekf->_H);
  }
  // Initialize the state transition function.
  ekf->f = options->f;
  // Initialize the measurement function.
  ekf->h = options->h;
  // Initialize the number of state variables.
  ekf->numberOfStates = options->n;
  // Initialize the use finite difference Jacobian flag.
  ekf->useFiniteDifferenceJacobian = options->useFiniteDifferenceJacobian;
  return EKF_SUCCESS;
} // EKFInit()

EKFReturnCodes EKFCleanup(EKFState *ekf)
{
  NULL_CHECK_EKF(ekf);

  if (ekf->mallocFlag)
  {
    FREE_MATRIX(ekf->_P);
    FREE_MATRIX(ekf->_K);
    FREE_MATRIX(ekf->_z);
    FREE_MATRIX(ekf->_F);
    FREE_MATRIX(ekf->_H);
    FREE_MATRIX(ekf->_F_TRANSPOSE);
    FREE_MATRIX(ekf->_H_TRANSPOSE);
  }
} // EKFCleanup()

// See EKF.h for documentation
EKFReturnCodes EKFPredict(EKFState *ekf)
{
  NULL_CHECK_EKF(ekf);
  // Use the system dynamics model to predict the state at the next time step.
  // Use the state transition function, f, to propagate the state: x_predicted = f(x).
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->A, ekf->x, ekf->x));
  // Calculate the Jacobian matrix, F, of the state transition function with respect to the state variables.
  if (ekf->useFiniteDifferenceJacobian)
  {
    calculateJacobian(ekf->x, ekf->_x_predicted, ekf->_F, ekf->f, ekf->numberOfStates);
  }
  // Update the state covariance matrix: P_predicted = _F * P * _F^T + Q.
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_F, ekf->P, ekf->P));
  MATRIX_MATH_RETURN_CHECK(tranposeMatrix(ekf->_F, ekf->_F_TRANSPOSE));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->P, ekf->_F_TRANSPOSE, ekf->P));
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->P, ekf->Q, ekf->P));
  return EKF_SUCCESS;
} // EKFPredict()

// See EKF.h for documentation
EKFReturnCodes EKFUpdate(EKFState *ekf, EKFMeasurement *measurement)
{
  NULL_CHECK_EKF(ekf);
  NULL_CHECK_EKF(measurement);
  // Calculate the Jacobian matrix, H, of the measurement function with respect to the state variables, evaluated at x_predicted.
  if (ekf->useFiniteDifferenceJacobian)
  {
    calculateJacobian(ekf->_x_predicted, ekf->_z, ekf->_H, ekf->h, ekf->numberOfStates);
  }
  // Calculate the Kalman gain: K = P_predicted * H^T * (H * P_predicted * H^T + R)^-1.
  MATRIX_MATH_RETURN_CHECK(tranposeMatrix(ekf->_H, ekf->_H_TRANSPOSE));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->P, ekf->_P));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->P, ekf->_H_TRANSPOSE, ekf->_K));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_H, ekf->P, ekf->P));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->P, ekf->_H_TRANSPOSE, ekf->P));
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->P, ekf->R, ekf->P));
  MATRIX_MATH_RETURN_CHECK(inverseMatrix(ekf->P, ekf->P));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->P, ekf->_K));
  // Calculate the measurement residual: y = z - h(x).
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_H, ekf->x, ekf->_z));
  MATRIX_MATH_RETURN_CHECK(subMatrix(measurement->z, ekf->_z, ekf->_z));
  // Update the state estimate: x = x + K * y.
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->_z, ekf->_z));
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->x, ekf->_z, ekf->x));
  // Update the state covariance matrix: P = (I - K * H) * P_predicted.
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->_H, ekf->P));
  MATRIX_MATH_RETURN_CHECK(idenityMatrixMinusA(ekf->P, ekf->P));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->P, ekf->_P, ekf->P));
} // EKFUpdate()

// ------------------------- Private Functions ------------------------- //
void calculateJacobian(EKFMatrix *x, EKFMatrix *_x_predicted, EKFMatrix *Jacobian, EKFStateTransitionFunction f, int numOfStates)
{
  int i, j;
  for (i = 0; i < numOfStates; ++i) {
    // Perturb the state variable by EPSILON.
    x->mat[i][0] += EPSILON;
    // Calculate the state transition function with the perturbed state variable.
    f(x, _x_predicted);
    // Calculate the Jacobian matrix, F, of the state transition function with respect to the state variables.
    for (j = 0; j < numOfStates; ++j) {
      Jacobian->mat[j][i] = (_x_predicted->mat[j][0] - x->mat[j][0]) / EPSILON;
    } // for (j = 0; j < numOfStates; ++j)
    // Reset the state variable.
    x->mat[i][0] -= EPSILON;
  } // for (i = 0; i < numOfStates; ++i
} // calculateJacobian()