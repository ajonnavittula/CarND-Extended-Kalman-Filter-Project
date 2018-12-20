#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj_(3,4);

  // Some variables and calculations for brevity
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  double C1 = px*px + py*py;
  double C2 = sqrt(C1);
  double C3 = vx*py - vy*px;
  double C4 = vy*px - vx*py;
  double C5 = pow(C2,3);

  // Jacobian calculation
  if(px == 0 && py == 0)
  {
    std::cout<<"Calculation error: Divide by zero";
  }
  else
  {
    Hj_ << px/C2, py/C2, 0, 0,
           -py/C1, px/C1, 0, 0,
           py*C3/C5, px*C4/C5, px/C2, py/C2;
  }

  return Hj_;
}
