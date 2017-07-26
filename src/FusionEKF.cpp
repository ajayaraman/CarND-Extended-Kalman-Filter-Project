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
              0, 0.0225;

  //measurement covariance matrix - radar
    
   //Radar's estimate of position is noisy compared to velocity and heading
   R_radar_ <<  1, 0, 0,
                0,  0.002, 0,
                0, 0,   0.1;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
             0, 1, 0, 0;

  //Pre - initialize F state transition matrix so we don't have to change all elements each time a packet is received
  MatrixXd F_in(4,4);
  F_in << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  VectorXd x_in(4);
  x_in << 0, 0, 0, 0;

  MatrixXd P_in = MatrixXd::Identity(4, 4);

  MatrixXd Q_in = MatrixXd::Identity(4, 4);

  ekf_.Init(x_in, P_in, F_in, H_laser_, R_laser_, Q_in);      
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  VectorXd meas = measurement_pack.raw_measurements_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //rho * cos(phi), rho * sin(phi)
      ekf_.x_(0) = meas(0) * cos(meas(1));
      ekf_.x_(1) = meas(0) * sin(meas(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // x, y 
      ekf_.x_(0) = meas(0);
      ekf_.x_(0) = meas(1);
    }

    //setup previous timestamp to be the timestamp of the first packet so that dt is initialized correctly
    previous_timestamp_ = measurement_pack.timestamp_;
   
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

  previous_timestamp_ = measurement_pack.timestamp_;
    
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
    

  //Compute Q process modeling noise
  float dt_2 = dt* dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  float noise_ax = 9;
  float noise_ay = 9;

  ekf_.Q_ << dt_4 * noise_ax/4, 0, dt_3 * noise_ax/2, 0,
	          0, dt_4 * noise_ay/4, 0, dt_3 * noise_ay/2,
	          dt_3 * noise_ax/2, 0, dt_2 * noise_ax, 0,
	          0,   dt_3 * noise_ay/2, 0, dt_2 * noise_ay;
    
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      
    // Setup Measurement and Measurement Noise matrices for Radar
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(meas);

  } else {
    // Laser updates
  
    // Setup Measurement and Measurement Noise matrices for Laser
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(meas);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
