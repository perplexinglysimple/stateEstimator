#ifndef Kalman_Filter
#define Kalman_Filter

//The matrices needed for the kalman filter
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

int kalmanFilterInit(struct kalmanMatrixes *);

int kalmanFilterPredict(struct kalmanMatrixes *);

int kalmanFilterUpdate(struct kalmanMatrixes *);

int kalmanFilterUpdateEKF(struct kalmanMatrixes *);

#endif