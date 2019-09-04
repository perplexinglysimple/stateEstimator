#include "matrixMath.h"


int kalmanFilterInit(struct *kalmanMatrixes) {
	
}

int kalmanFilterPredict(struct *kalmanMatrixes) {
	int error = multMatrix(kalmanMatrixes->x_, kalmanMatrixes->F_, kalmanMatrixes->TEMP_X_);
	if(error > 0) {
		return error;
	}
	error = tranposeMatrix(kalmanMatrixes->F_, kalmanMatrixes->F_TRANSPOSE);
	if(error > 0) {
		return error;
	}
	error = multMatrix(kalmanMatrixes->F_, kalmanMatrixes->P_, kalmanMatrixes->TEMP_FP_);
	if(error > 0) {
		return error;
	}
	error = multMatrix(kalmanMatrixes->TEMP_FP_,  kalmanMatrixes->F_TRANSPOSE, kalmanMatrixes->TEMP_FPFT_);
	if(error > 0) {
		return error;
	}
	error = addMatrix(kalmanMatrixes->TEMP_FPFT_, kalmanMatrixes->Q_, kalmanMatrixes->P_);
	if(error > 0) {
		return error;
	}
	return 0;
}

int kalmanFilterUpdate(struct *kalmanMatrixes) {
	//Transpose H_
	int error = ranposeMatrix(kalmanMatrixes->H_, kalmanMatrixes->H_TRANSPOSE);
	if(error > 0) {
		return error;
	}
	//PHT_
	error = multMatrix(kalmanMatrixes->P_, kalmanMatrixes->H_TRANSPOSE, kalmanMatrixes->TEMP_PHT_);
	if(error > 0) {
		return error;
	}
	
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  //Update State
  x_ = x_ + (K * y);
  //Update covariance matrix
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);  
  P_ = (I - K*H_) * P_;
}

int kalmanFilterUpdateEKF(struct *kalmanMatrixes) {
	
}

