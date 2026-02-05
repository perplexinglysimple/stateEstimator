#include "EKF.h"

// ------------------------- Private Function Prototypes ------------------------- //
/// @brief Calculate the Jacobian matrix of a function with respect to the state variables.
void calculateJacobian(EKFMatrix *x, EKFMatrix *_x_predicted, EKFMatrix *Jacobian, EKFStateTransitionFunction f, EKFState *ekf, EKFMatrix *baseline, void* userData);

// ------------------------- Public Functions ------------------------- //
// See EKF.h for documentation
EKFReturnCodes EKFInit(EKFState *ekf, EKFConfigOptions *options)
{
  LOG_FUNCTION();

  NULL_CHECK_EKF(ekf);
  NULL_CHECK_EKF(options);

  // Check to see if the options struct is initialized.
  NULL_CHECK_MATRIX(options->x0);
  NULL_CHECK_MATRIX(options->P0);
  NULL_CHECK_MATRIX(options->Q);
  NULL_CHECK_MATRIX(options->R);
  NULL_CHECK_MATRIX(options->A);

  // Check to see if the ekf struct is initialized.
  NULL_CHECK_MATRIX(ekf->x);
  NULL_CHECK_MATRIX(ekf->P);
  NULL_CHECK_MATRIX(ekf->Q);
  NULL_CHECK_MATRIX(ekf->R);
  NULL_CHECK_MATRIX(ekf->A);

  int numMeasurements = options->numberOfMeasurements;
  if (numMeasurements <= 0)
  {
    numMeasurements = options->n;
  }
  if (numMeasurements <= 0)
  {
    LOG_ERROR("EKFInit() failed because numberOfMeasurements and n are not set.");
    return EKF_ERROR;
  }
  options->numberOfMeasurements = numMeasurements;

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
  if (options->R->row != numMeasurements || options->R->col != numMeasurements)
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
    INIT_MATRIX(ekf->_K, ekf->P->row, numMeasurements);
    INIT_MATRIX(ekf->_z, numMeasurements, 1);
    INIT_MATRIX(ekf->_F, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_H, numMeasurements, ekf->P->col);
    INIT_MATRIX(ekf->_F_TRANSPOSE, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_H_TRANSPOSE, ekf->P->col, numMeasurements);
    INIT_MATRIX(ekf->_x_predicted, ekf->x->row, ekf->x->col);
    INIT_MATRIX(ekf->_S, numMeasurements, numMeasurements);
    INIT_MATRIX(ekf->_I, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_TEMP1, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_TEMP2, ekf->P->row, ekf->P->col);
    INIT_MATRIX(ekf->_TEMP3, ekf->x->row, ekf->x->col);
    INIT_MATRIX(ekf->_TEMP4, numMeasurements, ekf->P->col);
    INIT_MATRIX(ekf->_TEMP5, ekf->P->row, numMeasurements);
    INIT_MATRIX(ekf->_TEMP6, numMeasurements, numMeasurements);
    INIT_MATRIX(ekf->_TEMP7, numMeasurements, 1);
    INIT_MATRIX(ekf->_TEMP8, ekf->x->row, ekf->x->col);
  }
  else
  {
    // Check that the temporary storage matrices are the correct size and not null
    NULL_CHECK_MATRIX(ekf->_P);
    NULL_CHECK_MATRIX(ekf->_K);
    NULL_CHECK_MATRIX(ekf->_z);
    NULL_CHECK_MATRIX(ekf->_F);
    NULL_CHECK_MATRIX(ekf->_H);
    NULL_CHECK_MATRIX(ekf->_F_TRANSPOSE);
    NULL_CHECK_MATRIX(ekf->_H_TRANSPOSE);
    NULL_CHECK_MATRIX(ekf->_x_predicted);
    NULL_CHECK_MATRIX(ekf->_S);
    NULL_CHECK_MATRIX(ekf->_I);
    NULL_CHECK_MATRIX(ekf->_TEMP1);
    NULL_CHECK_MATRIX(ekf->_TEMP2);
    NULL_CHECK_MATRIX(ekf->_TEMP3);
    NULL_CHECK_MATRIX(ekf->_TEMP4);
    NULL_CHECK_MATRIX(ekf->_TEMP5);
    NULL_CHECK_MATRIX(ekf->_TEMP6);
    NULL_CHECK_MATRIX(ekf->_TEMP7);
    NULL_CHECK_MATRIX(ekf->_TEMP8);
  }
  // Initialize the state transition function.
  ekf->f = options->f;
  // Initialize the measurement function.
  ekf->h = options->h;
  // Initialize Jacobian callbacks (optional).
  ekf->jacobianF = options->jacobianF;
  ekf->jacobianH = options->jacobianH;
  // Initialize the state transition matrix function.
  ekf->updateAMatrix = options->updateAMatrix;
  // Initialize the number of state variables.
  ekf->numberOfStates = options->n;
  // Initialize the number of measurement variables.
  ekf->numberOfMeasurements = numMeasurements;
  // Initialize the use finite difference Jacobian flag.
  ekf->useFiniteDifferenceJacobian = options->useFiniteDifferenceJacobian;
  if (!ekf->useFiniteDifferenceJacobian && (!ekf->jacobianF || !ekf->jacobianH))
  {
    LOG_ERROR("EKFInit() requires jacobianF and jacobianH when useFiniteDifferenceJacobian is false.");
    return EKF_ERROR;
  }
  return EKF_SUCCESS;
} // EKFInit()

EKFReturnCodes EKFCleanup(EKFState *ekf)
{
  LOG_FUNCTION();

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
    FREE_MATRIX(ekf->_x_predicted);
    FREE_MATRIX(ekf->_S);
    FREE_MATRIX(ekf->_I);
    FREE_MATRIX(ekf->_TEMP1);
    FREE_MATRIX(ekf->_TEMP2);
    FREE_MATRIX(ekf->_TEMP3);
    FREE_MATRIX(ekf->_TEMP4);
    FREE_MATRIX(ekf->_TEMP5);
    FREE_MATRIX(ekf->_TEMP6);
    FREE_MATRIX(ekf->_TEMP7);
    FREE_MATRIX(ekf->_TEMP8);
  }

  return EKF_SUCCESS;
} // EKFCleanup()

// See EKF.h for documentation
EKFReturnCodes EKFPredict(EKFState *ekf, double time, void* userData)
{
  LOG_FUNCTION();

  NULL_CHECK_EKF(ekf);
  // Use the system dynamics model to predict the state at the next time step.
  // Get the new A vector.
  ekf->updateAMatrix(ekf->A, ekf->x, ekf, time, userData);
  // Propagate the state: x_predicted = f(x).
  ekf->f(ekf->x, ekf->_x_predicted, ekf, userData);
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_x_predicted, ekf->x));
  // Calculate the Jacobian matrix, F, of the state transition function with respect to the state variables.
  if (ekf->useFiniteDifferenceJacobian)
  {
    calculateJacobian(ekf->x, ekf->_x_predicted, ekf->_F, ekf->f, ekf, ekf->_TEMP3, userData);
  }
  else if (ekf->jacobianF)
  {
    ekf->jacobianF(ekf->x, ekf->_F, ekf, userData);
  }
  // Update the state covariance matrix: P_predicted = _F * P * _F^T + Q.
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_F, ekf->P, ekf->_TEMP1));
  MATRIX_MATH_RETURN_CHECK(transposeMatrix(ekf->_F, ekf->_F_TRANSPOSE));
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_TEMP1, ekf->_F_TRANSPOSE, ekf->P));
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->P, ekf->Q, ekf->_TEMP1));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_TEMP1, ekf->P));
  return EKF_SUCCESS;
} // EKFPredict()

// See EKF.h for documentation
EKFReturnCodes EKFUpdate(EKFState *ekf, EKFMeasurement *measurement)
{
  LOG_FUNCTION();

  NULL_CHECK_EKF(ekf);
  NULL_CHECK_EKF(measurement);
  if (measurement->z->row != ekf->numberOfMeasurements || measurement->z->col != 1)
  {
    LOG_ERROR("EKFUpdate() failed because measurement z is not the correct size.");
    return EKF_ERROR;
  }
  // Calculate the Jacobian matrix, H, of the measurement function with respect to the state variables, evaluated at _x_predicted.
  if (ekf->useFiniteDifferenceJacobian)
  {
    calculateJacobian(ekf->_x_predicted, ekf->_z, ekf->_H, ekf->h, ekf, ekf->_TEMP7, NULL);
  }
  else if (ekf->jacobianH)
  {
    ekf->jacobianH(ekf->_x_predicted, ekf->_H, ekf, NULL);
  }
  // Calculate the Kalman gain: K = P_predicted * H^T * (H * P_predicted * H^T + R)^-1.
  MATRIX_MATH_RETURN_CHECK(transposeMatrix(ekf->_H, ekf->_H_TRANSPOSE));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->P, ekf->_P)); // _P holds P_predicted

  // S = H * P_predicted * H^T + R
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_H, ekf->_P, ekf->_TEMP4)); // (m x n)*(n x n) = (m x n)
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_TEMP4, ekf->_H_TRANSPOSE, ekf->_S)); // (m x n)*(n x m) = (m x m)
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->_S, ekf->R, ekf->_TEMP6)); // (m x m)
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_TEMP6, ekf->_S));

  // Invert S with jitter if needed
  matrixReturnCodes invRet = inverseMatrixWithJitter(ekf->_S, ekf->_TEMP6, (matrixType)1e-6, 3, (matrixType)100);
  if (invRet != MATRIX_SUCCESS)
  {
    LOG_ERROR("EKFUpdate() failed to invert innovation matrix S.");
    return EKF_ERROR;
  }

  // K = P_predicted * H^T * S^-1
  MATRIX_MATH_RETURN_CHECK(transposeMatrix(ekf->_TEMP4, ekf->_TEMP5)); // (n x m) = (m x n)^T
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_TEMP5, ekf->_TEMP6, ekf->_K)); // (n x m)*(m x m) = (n x m)

  // Calculate the measurement residual: y = z - h(x_predicted).
  ekf->h(ekf->_x_predicted, ekf->_z, ekf, NULL);
  MATRIX_MATH_RETURN_CHECK(subMatrix(measurement->z, ekf->_z, ekf->_TEMP7));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_TEMP7, ekf->_z));

  // Update the state estimate: x = x + K * y.
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->_z, ekf->_TEMP3)); // (n x m)*(m x 1) = (n x 1)
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->x, ekf->_TEMP3, ekf->_TEMP8));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_TEMP8, ekf->x));

  // Joseph form covariance update: P = (I - K H) P (I - K H)^T + K R K^T
  MATRIX_MATH_RETURN_CHECK(setIdentityMatrix(ekf->_I));

  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->_H, ekf->_TEMP1)); // (n x m)*(m x n) = (n x n)
  MATRIX_MATH_RETURN_CHECK(subMatrix(ekf->_I, ekf->_TEMP1, ekf->_TEMP2)); // TEMP2 = I - K*H
  MATRIX_MATH_RETURN_CHECK(transposeMatrix(ekf->_TEMP2, ekf->_TEMP1)); // TEMP1 = (I - K*H)^T
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_TEMP2, ekf->P, ekf->_P)); // _P = (I-KH)*P_pred
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_P, ekf->_TEMP1, ekf->P)); // P = (I-KH)*P_pred*(I-KH)^T

  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_K, ekf->R, ekf->_TEMP5)); // (n x m)*(m x m) = (n x m)
  MATRIX_MATH_RETURN_CHECK(transposeMatrix(ekf->_K, ekf->_TEMP4)); // (m x n)
  MATRIX_MATH_RETURN_CHECK(multMatrix(ekf->_TEMP5, ekf->_TEMP4, ekf->_TEMP1)); // (n x m)*(m x n) = (n x n)
  MATRIX_MATH_RETURN_CHECK(addMatrix(ekf->P, ekf->_TEMP1, ekf->_TEMP2));
  MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf->_TEMP2, ekf->P));
  return EKF_SUCCESS;
} // EKFUpdate()

// ------------------------- Private Functions ------------------------- //
void calculateJacobian(EKFMatrix *x, EKFMatrix *_x_predicted, EKFMatrix *Jacobian, EKFStateTransitionFunction f, EKFState *ekf, EKFMatrix *baseline, void* userData)
{
  LOG_FUNCTION();

  int i, j;
  matrixType temp;
  matrixType epislon = EPSILON;
  int numOfStates = x->row;
  int numOutputs = Jacobian->row;
  // Compute baseline f(x) into provided baseline buffer.
  f(x, baseline, ekf, userData);
  for (i = 0; i < numOfStates; ++i) {
    // Perturb the state variable by EPSILON.
    SET_MATRIX(*x, i, 0, ACCESS_MATRIX(*x, i, 0) + epislon);
    // Calculate the state transition function with the perturbed state variable.
    // This takes in the state vector, ekf object, user data and outputs the predicted state vector.
    // Most of the time this will be a simple multiplication operation of the state transition matrix and the state vector.
    f(x, _x_predicted, ekf, userData);
    // Calculate the Jacobian matrix, F, of the state transition function with respect to the state variables.
    for (j = 0; j < numOutputs; ++j) {
      temp = (ACCESS_MATRIX(*_x_predicted, j, 0) - ACCESS_MATRIX(*baseline, j, 0)) / epislon;
      SET_MATRIX(*Jacobian, j, i, temp);
    } // for (j = 0; j < numOutputs; ++j)
    // Reset the state variable.
    SET_MATRIX(*x, i, 0, ACCESS_MATRIX(*x, i, 0) - epislon);
  } // for (i = 0; i < numOfStates; ++i
} // calculateJacobian()
