#include "EKF.h"

// System dynamics model function prototype
void TransistionFunction(EKFMatrix* x, EKFMatrix* x_predicted);
// Measurement function prototype
void MeasurementFunction(EKFMatrix* x, EKFMatrix* z);

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
 * A = [1 0 0.1  0    0    0   ]
 *     [0 1 0    0.1  0    0   ]
 *     [0 0 1    0    0.1  0   ]
 *     [0 0 0    1    0    0.1 ]
 *     [0 0 0    0    1    0   ]
 *     [0 0 0    0    0    1   ]
 * 
 * The tranistion function will be:
 * x_predicted = A*x
 * 
 * The measurement function will be:
 * z = x + noise and we cannot measure the velocity or acceleration
 */

int main()
{
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    EKFMeasurement z = {0};
    ekfType xadd, yadd = 0;
    int count = 0;
    bool continueFlag = true;
    // Initialize the options struct.s
    STATIC_MATRIX_DIRECTIVE(options.x0, 6, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, 6, 6, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, 6, 6, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, 6, 6, R);
    STATIC_MATRIX_DIRECTIVE(options.A, 6, 6, A);
    options.n = 6;
    options.f = TransistionFunction;
    options.h = MeasurementFunction;
    options.numberOfStates = 6;
    options.useFiniteDifferenceJacobian = true;
    options.mallocFlag = true;

    // Initialize the ekf struct.
    PRE_INIT_ALLOC(&ekf, 6, true);

    // Initialize the measurement struct.
    STATIC_MATRIX_DIRECTIVE(z.z, 6, 1, z);
    // This should succeed because the options struct has been initialized.
    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed to initialize the EKF.");
        return -1;
    }

    LOG_INFO("Start");

    // Run the EKF.
    while (continueFlag)
    {
        // Predict the next state.
        if (EKFPredict(&ekf) != EKF_SUCCESS)
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
        ekfType randNum = (ekfType)rand() / (ekfType)RAND_MAX;
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

void TransistionFunction(EKFMatrix *x, EKFMatrix *x_predicted)
{
    LOG_INFO("TransistionFunction() called.");
    struct matrix *A = NULL;
    STATIC_MATRIX_DIRECTIVE(A, 6, 6, dumdumvar);
    // Fill in the A matrix.
    ACCESS_STATIC_MATRIX(*A, 0, 1) = 1;
    ACCESS_STATIC_MATRIX(*A, 0, 2) = 0.1;
    ACCESS_STATIC_MATRIX(*A, 1, 1) = 1;
    ACCESS_STATIC_MATRIX(*A, 1, 3) = 0.1;
    ACCESS_STATIC_MATRIX(*A, 2, 2) = 1;
    ACCESS_STATIC_MATRIX(*A, 3, 3) = 1;
    ACCESS_STATIC_MATRIX(*A, 4, 4) = 1;
    ACCESS_STATIC_MATRIX(*A, 5, 5) = 1;
    // x_predicted = A*x
    multMatrix(A, x, x_predicted);
}

void MeasurementFunction(EKFMatrix *x, EKFMatrix *z)
{
    LOG_INFO("MeasurementFunction() called.");
    // The measurement will be the state plus some noise.
    // z = x
    copyMatrix(x, z);
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
