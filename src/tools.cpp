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
    * TODO:
    *   Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        // Ignore data if  it's not complete
        cout << " CalculateRMSE : Invalid estimations or  ground_truth vector size" << endl;
        return rmse;
    }
    
    for (int idx = 0; idx < estimations.size(); ++idx) {
        Eigen::VectorXd diff = estimations[idx] - ground_truth[idx];
        diff  = diff.array() * diff.array();
        rmse += diff;
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   *  Calculate a Jacobian here.
   */
    MatrixXd Hj(3,4);

    if (x_state.size() != 4) {
      cout << "CalculateJacobian () - Error - Statue vector has wrong size" << endl;
      return Hj;
    }

    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
    
    // Pre-compute values
    double term1 = px*px+py*py;
    double term2 = sqrt(term1);
    double term3 = (term1*term2);
    
    // check division by zero
    if (term1 < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    // compute the Jacobian matrix
    Hj <<  (px/term2),                (py/term2),               0,        0,
          -(py/term1),                (px/term1),               0,        0,
            py*(py*vx - px*vy)/term3, px*(px*vy - py*vx)/term3, px/term2, py/term2;

    return Hj;
}
