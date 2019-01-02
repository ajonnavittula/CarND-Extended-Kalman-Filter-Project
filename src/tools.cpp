#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //Used the RMSE code from EKF lectures

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()) {
    cout << "Size error: Unequal estimations or ground_truth" << endl;
    return rmse;
  }
  else if (estimations.size() == 0)
  {
    cout<< "Size error: Empty estimations vector"<<endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) 
  {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  MatrixXd Hj_(3,4);

  // Some variables and calculations for brevity
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float C1 = px*px + py*py;
  float C2 = sqrt(C1);
  float C3 = vx*py - vy*px;
  float C4 = vy*px - vx*py;
  float C5 = pow(C2,3);

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
