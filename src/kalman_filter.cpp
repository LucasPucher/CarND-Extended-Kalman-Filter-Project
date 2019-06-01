#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;			/* System state 					*/
	P_ = P_in;			/* Covariance matrix of system		*/
	F_ = F_in;			/* System model matrix				*/
	H_ = H_in;			/* Measurement matrix				*/
	R_ = R_in;			/* Measurement noise / covariance 	*/
	Q_ = Q_in;			/* Process noise / covariance		*/
}

void KalmanFilter::Predict() {
	x_ = F_ * x_;							/* Predict system state based on the system model matrix 			*/
	P_ = F_ * P_ * F_.transpose() + Q_;		/* Predict system state covariance based on the system model matrix */
}

void KalmanFilter::Update(const VectorXd &z) {
	
	/* Calculate some matrixes for speed optimization */
	MatrixXd Ht = H_.transpose();
	
	/* Calculation of Kalman gain	*/
	VectorXd y = z - H_ * x_;				/* Innovation, meaning difference between measurement and estimate 					*/
	MatrixXd S = H_ * P_ * Ht + R_;			/* Covariance of the innovation														*/
	MatrixXd K = P_ * Ht * S.inverse();		/* Calculation of the Kalman gain, meaning how much we can trust the measurement 	*/


	/* Refine estimate with measurement results */
	x_ = x_ + (K * y);						/* Update the new estimate with the measurement values 	*/	
	long i = x_.size();						/* Prepare creation of identity matrix 					*/
	MatrixXd I = MatrixXd::Identity(i, i);	/* Create identity matrix								*/
	P_ = (I - K * H_) * P_;					/* Calculate new covariance of the state 				*/
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/* Calculate some matrixes for speed optimization */
	MatrixXd Ht = H_.transpose();
	
	/* Convert prior estimate to radar format 										*/
	/* State equations are cartesian and radar (measurement equations) are polar 	*/
	/* We do not use a Jacobian here, but directly calculate the magnitudes			*/
	float range = sqrt(x_(0)*x_(0) + x_(1)*x_(1));	
	float bearing;
	if( (x_(0) == 0) && (x_(1) == 0))
	{
		bearing = 0;
		std::cout << "Undefined atan2!\n";
	}
	else
	{
		bearing = atan2(x_(1), x_(0));
	}
	
	float range_rate;
	/* Avoid division by zero	*/
	if (fabs(range) < 0.0001) {
	range_rate = 0;
	} 
	else {
	range_rate = (x_(0)*x_(2) + x_(1)*x_(3))/range;
	}
	VectorXd z_pred(3);
	z_pred << range, bearing, range_rate;
	
	/* Calculation of Kalman gain	*/
	VectorXd y = z - z_pred;
	
	/* Adjust for overflow			*/
	if (y(1) > M_PI) {
		y(1) -= 2 * M_PI;
	}
	else if(y(1) < (-M_PI)) {
		y(1) += 2 * M_PI;
	}
	
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd K = P_ * Ht * S.inverse();


	/* Refine estimate with measurement results */
	x_ = x_ + (K * y);						/* Update the new estimate with the measurement values 	*/	
	long i = x_.size();						/* Prepare creation of identity matrix 					*/
	MatrixXd I = MatrixXd::Identity(i, i);	/* Create identity matrix								*/
	P_ = (I - K * H_) * P_;					/* Calculate new covariance of the state 				*/
}	