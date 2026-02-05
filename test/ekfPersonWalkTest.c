#include "EKF.h"

// System dynamics model function prototype
void TransitionFunction(EKFMatrix* x, EKFMatrix* x_predicted, EKFState* ekf, void* userData);
// Measurement function prototype
void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
// Jacobian function prototypes
void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);

void MotionModelUpdateCallback(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);


int nanCheckVariable(EKFState *ekf);

/** Here is a picture of the state struct for a random walk model:
 *
 * x = [x]
 *     [y]
 *     [x_v]
 *     [y_v]
 *     [x_a]
 *     [y_a]
 * 
 * P0 is the initial state covariance matrix.
 * P0 = [1 0 0 0 0 0]
 *      [0 1 0 0 0 0]
 *      [0 0 1 0 0 0]
 *      [0 0 0 1 0 0]
 *      [0 0 0 0 1 0]
 *      [0 0 0 0 0 1]
 * 
 * Q is the process noise covariance matrix.
 * Q = [0.01 0    0    0    0    0   ]
 *     [0    0.01 0    0    0    0   ]
 *     [0    0    0.01 0    0    0   ]
 *     [0    0    0    0.01 0    0   ] 
 *     [0    0    0    0    0.01 0   ]
 *     [0    0    0    0    0    0.01]
 * 
 * R is the measurement noise covariance matrix.
 * R = [0.01 0    0    0    0    0   ]
 *     [0    0.01 0    0    0    0   ]
 *     [0    0    0.01 0    0    0   ]
 *     [0    0    0    0.01 0    0   ]
 *     [0    0    0    0    0.01 0   ]
 *     [0    0    0    0    0    0.01]
 * 
 * A is the state transition matrix.
 * A = [1 0 dt  0    0.5*dt^2    0   ]
 *     [0 1 0    dt  0    0.5*dt^2 ]
 *     [0 0 1    0    dt  0   ]
 *     [0 0 0    1    0    dt ]
 *     [0 0 0    0    1    0   ]
 *     [0 0 0    0    0    1   ]
 * 
 * The inputs to the ekf
 * 1. GPS position with LAT and LON and altitude. The LAT and LON are in the WGS84 coordinate system.
 *    x = [LAT]
 *        [LON]
 *        [ALT]
 *    GPS accuracy parameters
 *    DecayFactor = 0.5 (This is the random walk noise parameter)
 *    HorizontalPositionAccuracy = 1 (This is the horizontal accuracy of the GPS measurement)
 *    VerticalPositionAccuracy = 1 (This is the vertical accuracy of the GPS measurement)
 *    VelocityAccuracy = 1 (This is the velocity accuracy of the GPS measurement)
 * 2. The IMU measurements
 *   x = [x_a]
 *       [y_a]
 *       [z_a]
 *       [x_g]
 *       [y_g]
 *       [z_g]
 *       [x_m]
 *       [y_m]
 *       [z_m]
 *   IMU accuracy parameters
 *      Acclerometer
 *         MeasurmentRange = 2 (This is the measurement range of the accelerometer)
 *         Resolution = 0.1 (This is the resolution of the accelerometer)
 *         ConstantBias = 0.1 (This is the constant bias of the accelerometer)
 *         NoiseDensity = 0.1 (This is the noise density of the accelerometer)
 *     Gyroscope
 *        MeasurmentRange = 2 (This is the measurement range of the gyroscope)
 *        Resolution = 0.1 (This is the resolution of the gyroscope)
 *        ConstantBias = 0.1 (This is the constant bias of the gyroscope)
 *        AxesMisalignment = 0.1 (This is the axes misalignment of the gyroscope)
 *        NoiseDensity = 0.1 (This is the noise density of the gyroscope)
 *   Magnetometer
 *     MeasurmentRange = 2 (This is the measurement range of the magnetometer)
 *     Resolution = 0.1 (This is the resolution of the magnetometer)
 *     ConstantBias = 0.1 (This is the constant bias of the magnetometer)
 *     NoiseDensity = 0.1 (This is the noise density of the magnetometer)
 * 
 * This will result in the H matrix being:
 * 
 */

#define stateSize 6

int main()
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    EKFMeasurement z = {0};
    ekfType xadd, yadd = 0;
    int count = 0;
    bool continueFlag = true;
    int iterations = 0;
    const int maxIterations = 500;
    // Initialize the options struct.s
    STATIC_MATRIX_DIRECTIVE(options.x0, stateSize, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, stateSize, stateSize, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, stateSize, stateSize, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, stateSize, stateSize, R);
    STATIC_MATRIX_DIRECTIVE(options.A, stateSize, stateSize, A);
    options.n = 6;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.updateAMatrix = MotionModelUpdateCallback;
    options.numberOfStates = 6;
    options.numberOfMeasurements = 6;
    options.useFiniteDifferenceJacobian = false;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;
    options.mallocFlag = true;

    ekfType startx0[stateSize][1] = {{0}, {0}, {0}, {0}, {0}, {0}};
    ekfType startP0[stateSize][stateSize] = {{1, 0, 0, 0, 0, 0},
                                             {0, 1, 0, 0, 0, 0},
                                             {0, 0, 1, 0, 0, 0},
                                             {0, 0, 0, 1, 0, 0},
                                             {0, 0, 0, 0, 1, 0},
                                             {0, 0, 0, 0, 0, 1}};
    ekfType startQ[stateSize][stateSize] = {{0.01, 0, 0, 0, 0, 0},
                                            {0, 0.01, 0, 0, 0, 0},
                                            {0, 0, 0.01, 0, 0, 0},
                                            {0, 0, 0, 0.01, 0, 0},
                                            {0, 0, 0, 0, 0.01, 0},
                                            {0, 0, 0, 0, 0, 0.01}};
    ekfType startR[stateSize][stateSize] = {{0.01, 0, 0, 0, 0, 0},
                                            {0, 0.01, 0, 0, 0, 0},
                                            {0, 0, 0.01, 0, 0, 0},
                                            {0, 0, 0, 0.01, 0, 0},
                                            {0, 0, 0, 0, 0.01, 0},
                                            {0, 0, 0, 0, 0, 0.01}};
    ekfType startA[stateSize][stateSize] = {{1, 0, 0.1, 0, 0.005, 0},
                                            {0, 1, 0, 0.1, 0, 0.005},
                                            {0, 0, 1, 0, 0.1, 0},
                                            {0, 0, 0, 1, 0, 0.1},
                                            {0, 0, 0, 0, 1, 0},
                                            {0, 0, 0, 0, 0, 1}};
                                    
    COPY_2DARRAY_TO_MATRIX(startx0, options.x0);
    COPY_2DARRAY_TO_MATRIX(startP0, options.P0);
    COPY_2DARRAY_TO_MATRIX(startQ, options.Q);
    COPY_2DARRAY_TO_MATRIX(startR, options.R);
    COPY_2DARRAY_TO_MATRIX(startA, options.A);

    // Initialize the ekf struct.
    PRE_INIT_ALLOC(&ekf, 6, 6, true);

    // Initialize the measurement struct.
    STATIC_MATRIX_DIRECTIVE(z.z, 6, 1, z);
    // This should succeed because the options struct has been initialized.
    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to initialize the EKF.");
        return -1;
    }

    LOG_FUNCTION();

    // Run the EKF.
    while (continueFlag)
    {
        if (iterations++ >= maxIterations)
        {
            LOG_INFO("End");
            LOG_INFO("Reached max iterations without leaving bounds.");
            break;
        }
        // TODO need to add "action" into predict function so that we can update the non-linear motion model with it (A)
        // Predict the next state.
        if (EKFPredict(&ekf, .1, NULL) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFPredict() failed to predict the next state.");
            return -1;
        }
        if (nanCheckVariable(&ekf) != 0)
        {
            LOG_ERROR("nanCheckVariable() found a NaN in the EKFState struct.");
            return -1;
        }
        LOG_INFO("Predicted state:");
        LOG_INFO("x = %f, y = %f, x_v = %f, y_v = %f, x_a = %f, y_a = %f", ACCESS_MATRIX(*(ekf.x), 0, 0), ACCESS_MATRIX(*(ekf.x), 1, 0), ACCESS_MATRIX(*(ekf.x), 2, 0), ACCESS_MATRIX(*(ekf.x), 3, 0), ACCESS_MATRIX(*(ekf.x), 4, 0), ACCESS_MATRIX(*(ekf.x), 5, 0));
        // Only update the state with a measurement every 10 iterations.
        if (count % 10 != 0)
        {
            count += 1;
            continue;
        }
        // Generate a measurement.
        // There is a 10% change that the measurement will have a large error.
        // This is to simulate a sensor failure.
        MATRIX_MATH_RETURN_CHECK(copyMatrix(ekf.x, z.z));
        ACCESS_STATIC_MATRIX(*(z.z), 2, 0) = 0;
        ACCESS_STATIC_MATRIX(*(z.z), 3, 0) = 0;
        ACCESS_STATIC_MATRIX(*(z.z), 4, 0) = 0;
        ACCESS_STATIC_MATRIX(*(z.z), 5, 0) = 0;
        // Generate a random number between 0 and 1.
        ekfType randNum = .9;//(ekfType)rand() / (ekfType)RAND_MAX;
        // If the random number is less than 0.1, then the measurement will have a large error.
        if (randNum < 0.1)
        {
            // We can only move by 0.1 in any direction.
            xadd = (ekfType)rand() / (ekfType)RAND_MAX * 10;
            yadd = (ekfType)rand() / (ekfType)RAND_MAX * 10;
        }
        else
        {
            // We can only move by 0.1 in any direction.
            xadd = (ekfType)rand() / (ekfType)RAND_MAX * 0.1;
            yadd = (ekfType)rand() / (ekfType)RAND_MAX * 0.1;
        }
        ACCESS_STATIC_MATRIX(*(z.z), 0, 0) += xadd;
        ACCESS_STATIC_MATRIX(*(z.z), 1, 0) += yadd;
        LOG_INFO("Measurement:");
        LOG_INFO("x = %f, y = %f, x_v = %f, y_v = %f, x_a = %f, y_a = %f", ACCESS_MATRIX(*(z.z), 0, 0), ACCESS_MATRIX(*(z.z), 1, 0), ACCESS_MATRIX(*(z.z), 2, 0), ACCESS_MATRIX(*(z.z), 3, 0), ACCESS_MATRIX(*(z.z), 4, 0), ACCESS_MATRIX(*(z.z), 5, 0));
        // Update the state with a measurement.
        if (EKFUpdate(&ekf, &z) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFUpdate() failed to update the state with a measurement.");
            return -1;
        }
        LOG_INFO("Updated state:");
        LOG_INFO("x = %f, y = %f, x_v = %f, y_v = %f, x_a = %f, y_a = %f", ACCESS_MATRIX(*(ekf.x), 0, 0), ACCESS_MATRIX(*(ekf.x), 1, 0), ACCESS_MATRIX(*(ekf.x), 2, 0), ACCESS_MATRIX(*(ekf.x), 3, 0), ACCESS_MATRIX(*(ekf.x), 4, 0), ACCESS_MATRIX(*(ekf.x), 5, 0));
        if (nanCheckVariable(&ekf) != 0)
        {
            LOG_ERROR("nanCheckVariable() found a NaN in the EKFState struct.");
            return -1;
        }
        // Check to see if we should continue.
        if (ekf.x->mat[0][0] > 100 || ekf.x->mat[0][0] < -100 || ekf.x->mat[1][0] > 100 || ekf.x->mat[1][0] < -100)
        {
            continueFlag = false;
            LOG_INFO("End");
            LOG_INFO("The state has gone out of bounds.");
        }
        count = 0;
    }

    // Cleanup the EKF.
    if (EKFCleanup(&ekf) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFCleanup() failed to cleanup the EKF.");
        return -1;
    }
    return 0;
}

/**
 * This is proving too hard to get working and testing. Forcing the user to provide the Jacobian is the best option.
 * TODO: Add the jacobi funtion to the struct and use it in the predict and update functions.
 * where it would be like this
 * r = sqrt (x(1)^2+x(3)^2);
 * b = atan2(x(3),x(1));
 * H = [ cos(b) 0 sin(b) 0;
        -sin(b)/r 0 cos(b)/r 0];]
 */
void TransitionFunction(EKFMatrix *x, EKFMatrix *x_predicted, EKFState* ekf, void* userData)
{
    (void)userData;
    LOG_FUNCTION();
    
    // x_predicted = A*x. Using A from the EKFState struct.
    multMatrix(ekf->A, x, x_predicted);
}

void MeasurementFunction(EKFMatrix *x, EKFMatrix *z, EKFState* ekf, void* userData)
{
    (void)ekf;
    (void)userData;
    LOG_FUNCTION();
    // The measu
    copyMatrix(x, z);
}

void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)userData;
    LOG_FUNCTION();

    for (int i = 0; i < J->row; ++i)
    {
        for (int j = 0; j < J->col; ++j)
        {
            if (J->jaggedAlloc)
            {
                J->mat[i][j] = ACCESS_MATRIX(*(ekf->A), i, j);
            }
            else
            {
                ACCESS_STATIC_MATRIX(*J, i, j) = ACCESS_MATRIX(*(ekf->A), i, j);
            }
        }
    }
}

void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void)x;
    (void)ekf;
    (void)userData;
    LOG_FUNCTION();

    for (int i = 0; i < J->row; ++i)
    {
        for (int j = 0; j < J->col; ++j)
        {
            matrixType val = (i == j) ? 1 : 0;
            if (J->jaggedAlloc)
            {
                J->mat[i][j] = val;
            }
            else
            {
                ACCESS_STATIC_MATRIX(*J, i, j) = val;
            }
        }
    }
}

void MotionModelUpdateCallback(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double timeElasped, void* userData)
{
    (void)x;
    (void)ekf;
    (void)userData;
    LOG_FUNCTION();
    // The A matrix is the same as the A matrix in the TransitionFunction.
    // Fill in the A matrix with the delta time and state.
    ekfType newA[stateSize][stateSize] = {{1, 0, timeElasped, 0, 0.5 * timeElasped * timeElasped, 0},
                                          {0, 1, 0, timeElasped, 0, 0.5 * timeElasped * timeElasped},
                                          {0, 0, 1, 0, timeElasped, 0},
                                          {0, 0, 0, 1, 0, timeElasped},
                                          {0, 0, 0, 0, 1, 0},
                                          {0, 0, 0, 0, 0, 1}};
    COPY_2DARRAY_TO_MATRIX(newA, A);
}

// We want to check all the variables in the EKFState struct to make sure they are not NaN.
int nanCheckVariable(EKFState *ekf)
{
    // Check the x matrix.
    if (nanCheckMatrix(ekf->x) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the x matrix.");
        return -1;
    }
    // Check the P matrix.
    if (nanCheckMatrix(ekf->P) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the P matrix.");
        return -1;
    }
    // Check the Q matrix.
    if (nanCheckMatrix(ekf->Q) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the Q matrix.");
        return -1;
    }

    // Check the R matrix.
    if (nanCheckMatrix(ekf->R) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the R matrix.");
        return -1;
    }

    // Check the A matrix.
    if (nanCheckMatrix(ekf->A) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the A matrix.");
        return -1;
    }

    // Check the _P matrix.
    if (nanCheckMatrix(ekf->_P) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the _P matrix.");
        return -1;
    }

    // Check the _K matrix.
    if (nanCheckMatrix(ekf->_K) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the _K matrix.");
        return -1;
    }

    // Check the _z matrix.
    if (nanCheckMatrix(ekf->_z) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the _z matrix.");
        return -1;
    }

    // Check the _F matrix.
    if (nanCheckMatrix(ekf->_F) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the _F matrix.");
        return -1;
    }

    // Check the _H matrix.
    if (nanCheckMatrix(ekf->_H) != 0)
    {
        LOG_ERROR("nanCheckVariable() found a NaN in the _H matrix.");
    }

    return 0;
}
