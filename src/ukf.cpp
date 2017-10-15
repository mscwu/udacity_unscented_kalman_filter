#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  // Number of states
  n_x_ = 5;
  
  // Number of augmented states
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  // Vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initial state covariance matrix
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  // Time stamp
  time_us_ = 0.0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/  	
  if (!is_initialized_) {
  	cout << "UKF: " << endl;

    // begin debug
    cout << "Initializing..." << endl;
    // end debug

  	x_ << 1, 1, 1, 1, 0.1;

  	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
  		/**
  		 * Convert radar from polar to cartesian coordinates and initiate state
  		 */
  		double rho = meas_package.raw_measurements_[0];
  		double phi = meas_package.raw_measurements_[1];
  		x_[0] = rho * cos(phi);
  		x_[1] = rho * sin(phi);
  	}
  	else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
  		/**
  		 * Initiate state with lidar data directly
  		 */
  		double px = meas_package.raw_measurements_[0];
  		double py = meas_package.raw_measurements_[1];
  		x_[0] = px;
  		x_[1] = py;
  	}
  	time_us_ =  meas_package.timestamp_;
  	// done initialization
  	is_initialized_ = true;
  	return;
  }
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // begin debug
  cout << "Predicting..." << endl;
  // end debug
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
  *  Update
  ****************************************************************************/
  // begin debug
  cout << "Updating..." << endl;
  // end debug

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // begin debug
    cout << "Received radar data..." << endl;
    // end debug

  	UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // begin debug
    cout << "Received lidar data..." << endl;
    // end debug

  	UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Set weights
  weights_.fill(0.0);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (unsigned int i = 1; i<2*n_aug_+1; ++i) {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
  }

  // Create augmented sigma points
  MatrixXd Xsig_aug = AugmentedSigmaPoints();

  // Create predicted sigma points
  SigmaPointsPrediction(delta_t, Xsig_aug);

  // Predict state mean
  VectorXd wt = weights_.transpose();
  x_.fill(0.0);
  x_ = Xsig_pred_ * wt; 
  
  // Predict state covariance matrix
  P_.fill(0.0);
  for (unsigned int i=0;i<2*n_aug_+1; ++i){
      VectorXd diff = Xsig_pred_.col(i) - x_;
      // normalize angle
      while(diff(3) > 2*M_PI) diff(3) -= 2 * M_PI;
      while(diff(3) < -2*M_PI) diff(3) += 2 * M_PI;
      P_ = P_ + weights_(i) * diff * (diff.transpose());
  }     
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // ground truth	
  VectorXd z = meas_package.raw_measurements_;

  // measurement size
  int n_z = z.size();

  // predicted sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // predicted measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);

  // begin debug
  cout << "Predicting lidar measurement..." << endl;
  // end debug

  PredictLidarMeasurement(&z_pred, &S, &Zsig);

  // begin debug
  cout << "Updating state and state covariance..." << endl;
  // end debug
   
  UpdateStateCovariance(n_z, Zsig, z_pred, S, z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // ground truth 
  VectorXd z = meas_package.raw_measurements_;

  // measurement size
  int n_z = z.size();

  // predicted sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // predicted measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);

  // begin debug
  cout << "Predicting radar measurement..." << endl;
  // end debug

  PredictRadarMeasurement(&z_pred, &S, &Zsig);

  // begin debug
  cout << "Updating state and state covariance..." << endl;
  // end debug

  UpdateStateCovariance(n_z, Zsig, z_pred, S, z);
}

MatrixXd UKF::AugmentedSigmaPoints() {
  // Create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Create augmented mean state
  x_aug << x_, 0, 0;

  // Create augmented covariance matrix
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0,
       0, std_yawdd_ * std_yawdd_;
  P_aug << P_, MatrixXd::Zero(n_x_, 2),
           MatrixXd::Zero(2, n_x_), Q; 

  // Create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  double root = sqrt(lambda_ + n_aug_);
  MatrixXd sqroot = root * A;
  MatrixXd neg_sqroot = -sqroot;

  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = (sqroot.colwise() + x_aug);
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) = (neg_sqroot.colwise() + x_aug); 

  // Return augmented sigma points
  return Xsig_aug;	
}

void UKF::SigmaPointsPrediction(double delta_t, const MatrixXd &Xsig_aug) {
  // Construct noise vector
  MatrixXd noise = MatrixXd(n_x_, 1);
  for (unsigned int i=0; i<Xsig_pred_.cols();++i) {
      // get states
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double phi = Xsig_aug(3,i);
      double phi_d = Xsig_aug(4,i);
      double nu_a = Xsig_aug(5,i);
      double nu_phidd = Xsig_aug(6,i);
      
      if (phi_d != 0) {
          Xsig_pred_(0, i) = px + v * (sin(phi + phi_d * delta_t) - sin(phi)) / phi_d;
          Xsig_pred_(1, i) = py + v * (-cos(phi + phi_d * delta_t) + cos(phi)) / phi_d;

      }
      else {
          Xsig_pred_(0, i) = px + v * cos(phi) * delta_t;
          Xsig_pred_(1, i) = py + v * sin(phi) * delta_t;
      }
        Xsig_pred_(2, i) = v;
        Xsig_pred_(3, i) = phi + phi_d * delta_t;
        Xsig_pred_(4, i) = phi_d;
        
        noise(0, 0) = 0.5 * delta_t * delta_t * cos(phi) * nu_a;
        noise(1, 0) = 0.5 * delta_t * delta_t * sin(phi) * nu_a;
        noise(2, 0) = delta_t * nu_a;
        noise(3, 0) = 0.5 * delta_t * delta_t * nu_phidd;
        noise(4, 0) = delta_t * nu_phidd;
        
        Xsig_pred_.col(i) += noise;
    }
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  // Transform sigma points into measurement space
  // Using radar measurement model

  // Create matrix for sigma points in measurement space
  const int n_z = 3;	
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);	

  for (unsigned int i=0;i<2*n_aug_+1;++i) {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      double rho = sqrt(px * px + py * py);
      double phi = atan2(py, px);
      double phi_d = (px * cos(psi) * v + py * sin(psi) * v) / rho;
      Zsig.col(i) <<rho, phi, phi_d;
  }

  // Predicted measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  z_pred = Zsig * weights_; 

  // Predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  S << PredictMeasurementCommon(n_z, Zsig, z_pred);

  // Measurement noise R
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;

  // Add noise to S
  S += R;

  // Write result
  *z_out = z_pred;
  *S_out = S;	
  *Zsig_out = Zsig;
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  // Transform sigma points into measurement space
  // Using lidar measurement model

  // Create matrix for sigma points in measurment space
  const int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  for (unsigned int i=0;i<2*n_aug_+1;++i) {
  	  double px = Xsig_pred_(0, i);
  	  double py = Xsig_pred_(1, i);
  	  Zsig.col(i) << px, py;
  }

  // Predicted measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  z_pred = Zsig * weights_;

  // Predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  S << PredictMeasurementCommon(n_z, Zsig, z_pred);

  // Measurement noise R
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  R(0,0) = std_laspx_ * std_laspx_;
  R(1,1) = std_laspy_ * std_laspy_;

  // Add noise to S
  S += R;

  // Write result
  *z_out = z_pred;
  *S_out = S;	
  *Zsig_out = Zsig;  
}

MatrixXd UKF::PredictMeasurementCommon(const int n_z, 
	                                     const MatrixXd &Zsig,
                                       const VectorXd &z_pred) {


  // Predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (unsigned int i=0;i<2*n_aug_+1;++i) {
    VectorXd zdiff = Zsig.col(i) - z_pred;

    // normalize angle
    while(zdiff(1)>M_PI) zdiff(1) -= 2. * M_PI;
    while(zdiff(1)<-M_PI) zdiff(1) += 2. * M_PI;

    S += weights_(i) * zdiff * (zdiff.transpose());
  }
  
  return S;	
}

void UKF::UpdateStateCovariance(const int n_z,
	                              const MatrixXd &Zsig,
	                              const VectorXd &z_pred,
	                              const MatrixXd &S,
	                              const VectorXd &z) {
  // Calculate cross correlation matrix

  // begin debug
  cout << "Computing cross correlation matrix..." << endl;
  // end debug

  MatrixXd Tc = MatrixXd(n_x_, n_z);	
  Tc.fill(0.0);
  for (unsigned int i=0;i<2*n_aug_+1;++i) {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      // normalize angle
      while(xdiff(3)>M_PI) xdiff(3) -= 2. * M_PI;
      while(xdiff(3)<-M_PI) xdiff(3) += 2. * M_PI;
      VectorXd zdiff = Zsig.col(i) - z_pred;
      while(zdiff(1)>M_PI) zdiff(1) -= 2. * M_PI;
      while(zdiff(1)<-M_PI) zdiff(1) += 2. * M_PI;      
      Tc += weights_(i) * xdiff * zdiff.transpose();
  }
  // calculate Kalman gain K;

  // begin debug
  cout << Tc << endl;
  cout << "Computing kalman gain K..." << endl;
  // end debug

  MatrixXd K = Tc * S.inverse();
  VectorXd zdiff = z - z_pred;
  while(zdiff(1)>M_PI) zdiff(1) -= 2. * M_PI;
  while(zdiff(1)<-M_PI) zdiff(1) += 2. * M_PI;

  // begin debug
  cout << "Update state mean and covariance matrix..." << endl;
  // end debug   
  // update state mean and covariance matrix
  x_ += K * zdiff;
  P_ -= K * S * K.transpose();	
}
