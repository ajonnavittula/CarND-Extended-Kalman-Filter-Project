#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
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
  
  // State prediction using KF equations
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Calculate Kalman gain

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //Calculate new estimated positions
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  // Calcualte EKF gain
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  if(px == 0 && py == 0)
  {
    cout<<"Calculation error: Divide by zero\n";
    return;
  }

  double h1 = sqrt(px*px + py*py);
  double h2 = atan2(py,px);
  double h3 = (px*vx + py*vy)/h1;
  VectorXd z_pred(3);
  z_pred << h1,
            h2,
            h3;
  cout<<"Created variables for EKF\n";
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();

  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;

  MatrixXd K = PHt * Si;
  cout<<"Estimated EKF gain\n";
  //New estimated positions
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  P_ = (I - K * H_) * P_;
  cout<<"Calcualted new states\n";
}
