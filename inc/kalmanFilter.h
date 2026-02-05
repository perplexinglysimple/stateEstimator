/**
 * @file kalmanFilter.h
 * @brief Linear Kalman filter state container and API.
 */
#ifndef Kalman_Filter
#define Kalman_Filter

/** @brief Matrices needed for the linear Kalman filter. */
struct kalmanMatrixes {
	struct matrix *x_;
	struct matrix *P_;
	struct matrix *F_;
	struct matrix *H_;
	struct matrix *R_;
	struct matrix *Q_;
	struct matrix *F_TRANSPOSE;
	struct matrix *H_TRANSPOSE;
	struct matrix *TEMP_X_;
	struct matrix *TEMP_FP_;
	struct matrix *TEMP_FPFT_;
	struct matrix *TEMP_PHT_;
};

/**
 * @brief Initialize a Kalman filter state.
 * @return 0 on success; non-zero on failure.
 */
int kalmanFilterInit(struct kalmanMatrixes *);

/**
 * @brief Predict step for the linear Kalman filter.
 * @return 0 on success; non-zero on failure.
 */
int kalmanFilterPredict(struct kalmanMatrixes *);

/**
 * @brief Update step for the linear Kalman filter.
 * @return 0 on success; non-zero on failure.
 */
int kalmanFilterUpdate(struct kalmanMatrixes *);

#endif
