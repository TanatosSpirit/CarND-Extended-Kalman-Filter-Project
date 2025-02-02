#include "kalman_filter.h"

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
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K =  P_ * H_.transpose() * S.inverse();

  MatrixXd I; // Identity matrix
  I = MatrixXd::Identity(4, 4);
  // new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  VectorXd x_polar = tools.Cart2Polar(x_);

  VectorXd y = z - x_polar;

  // Normalize phi to [-pi, pi]
  while (y(1) > M_PI)  y(1) -= 2.0 * M_PI;
  while (y(1) < -M_PI) y(1) += 2.0 * M_PI;

  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K =  P_ * H_.transpose() * S.inverse();

  MatrixXd I; // Identity matrix
  I = MatrixXd::Identity(4, 4);
  // new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
