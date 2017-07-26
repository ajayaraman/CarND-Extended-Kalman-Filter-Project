#include "kalman_filter.h"
#include "tools.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;

#define SMALL_NUMBER (0.000001F)

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  	// KF Measurement update step
  MatrixXd Ht = H_.transpose();
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
		
  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //Hj computed in FusionEKF and H is set to Hj in FusionEKF.

  //Project state to measurement space using EKF nonlinear h(x)
  VectorXd hx(3);
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float phi = atan2(py, px);
  
  //constrain to -pi to pi for KF
  while (phi < - PI)
  {
    phi = phi + 2 * PI;
  }

  while(phi > PI)
  {
    phi = phi - 2 * PI;
  }
  
  float rho = sqrt(px * px + py * py) + SMALL_NUMBER;

  float rho_dot = (px * vx + py* vy) / rho;
  //project
  hx << rho, phi, rho_dot;

  // error
  VectorXd y = z - hx;

  //update Kalman gain
  MatrixXd Hjt = H_.transpose();
  MatrixXd S = H_ * P_ * Hjt + R_;
  MatrixXd K = P_ * Hjt * S.inverse();

  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
