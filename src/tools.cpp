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
  // estimations size should be the same as ground_truth
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;
  const int n = estimations.size();
  if (n == 0 || n != ground_truth.size()) {
    cout << "Estimations must be non-empty and the same length as ground truth!" << endl;
  }
  else {
    for (unsigned int i=0; i<n; ++i) {
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array() * residual.array();
      rmse += residual;
    }

    rmse = rmse.array() / n;
    rmse = rmse.array().sqrt();    
  }
  return rmse;
    
}