#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float px2 = px * px;
  float py2 = py * py;
  float sum_px2py2 = px2 + py2;
  float sqrt_px2py2 = sqrt(sum_px2py2);
  float px2py2_1p5 = pow(sum_px2py2, 1.5);
    
  bool isDenZero = fabs(sum_px2py2) < 0.0001;
  //check division by zero
	if(isDenZero)
  {
    std::cout << "Cannot Calculate Jacobian";
    return Hj;
  }
	//compute the Jacobian matrix
  Hj << px/sqrt_px2py2, py/sqrt_px2py2, 0 , 0,
          -py/sum_px2py2, px/sum_px2py2, 0 , 0,
           py * (vx * py - vy * px)/px2py2_1p5, px * (vy * px - vx * py)/px2py2_1p5, px/sqrt_px2py2, py/sqrt_px2py2;
  return Hj;
}
