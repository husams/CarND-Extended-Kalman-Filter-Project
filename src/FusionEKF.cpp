#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  /**
    * TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
    */
    // Initializing P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0,    0,
               0, 1, 0,    0,
               0, 0, 1000, 0,
               0, 0, 0,    1000;
    
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
     * TODO:
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * Remember: you'll need to convert radar from polar to cartesian coordinates.
     */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        cout << "EKF : First measurement RADAR" << endl;
      /**
       * Convert radar from polar to cartesian coordinates and initialize state.
       */
        double rho     = measurement_pack.raw_measurements_[0]; // range
        double phi     = measurement_pack.raw_measurements_[1];
        double rho_dot = measurement_pack.raw_measurements_[2];
        
        double x = rho * cos(phi);
        double y = rho * sin(phi);
        double vx = rho_dot * cos(phi);
        double vy = rho_dot * sin(phi);
        
        /*
        x = (x < 0.0001)? 0.0001 : x;
        y = (y < 0.0001)? 0.0001 : y;
        */
        
        ekf_.x_ << x, y, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        cout << "EKF : First measurement LASER" << endl;
      /**
       * Initialize state.
       */
        double x = measurement_pack.raw_measurements_[0];
        double y = measurement_pack.raw_measurements_[1];
        ekf_.x_ << x, y, 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_ ;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
    * TODO:
    * Update the state transition matrix F according to the new elapsed time.
    *  - Time is measured in seconds.
    * Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    // save last timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // State transition matrix update
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
               0, 1, 0,  dt,
               0, 0, 1,  0,
               0, 0, 0,  1;
    
    // Noise covariance matrix computation
    double noise_ax = 9.0;
    double noise_ay = 9.0;

    double dt2 = pow(dt, 2);
    double c1  = pow(dt, 4) / 4.0;
    double c2  = pow(dt, 3) / 2.0;
    
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << c1*noise_ax,   0,            c2*noise_ax,   0,
               0,             c1*noise_ay,  0,             c2 * noise_ay,
               c2 * noise_ax, 0,            dt2*noise_ax,  0,
               0,             c2*noise_ay,  0,             dt2*noise_ay;
    
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   * TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
      // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
