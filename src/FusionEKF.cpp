#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = MatrixXd(3, 4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
              
  //Initialization of EKF variables
  ekf_.x_ = VectorXd(4);

  ekf_.P_ = MatrixXd(4,4);
  

  ekf_.F_ = MatrixXd(4,4);  


  ekf_.Q_ = MatrixXd(4,4);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */

  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      ekf_.x_ << rho*cos(phi),
                 rho*sin(phi),
                 rho_dot*cos(phi),
                 rho_dot*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      // Initialize position vector to measured values
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
    }

    // Set timestamp value to measurement timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    // done initializing, no need to predict or update

    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  // State transition matrix updation

  //Delta T in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;


  // Give noise variance values noise ax = 9 and noise ay = 9
  float noise_ax = 9;
  float noise_ay = 9;

  // Update State transition matrix with current dt
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // Noise covariance matrix updation
  ekf_.Q_ << (dt4*noise_ax)/4, 0, (dt3*noise_ax)/2, 0,
              0, (dt4*noise_ay)/4, 0, (dt3*noise_ay)/2,
              (dt3*noise_ax)/2, 0, dt2*noise_ax, 0,
              0, dt3*noise_ay/2, 0, dt2*noise_ay;
  //std::cout<<"Q Matrix completed\n";


  ekf_.Predict();
  //std::cout<<"EKF predicted\n";
  /**
   * Update
   */


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    //std::cout<<"Calculating Jacobian\n";
    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    //std::cout<<"Updating EKF for radar\n";
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    //std::cout<<"EKF Update complete\n";

  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    //std::cout<<"Updating KF\n";
    ekf_.Update(measurement_pack.raw_measurements_);
    //std::cout<<"KF updated\n";
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
