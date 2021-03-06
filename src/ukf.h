#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* NIS check for radar
  double e_radar_;

  ///* NIS check for lidar
  double e_laser_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * Create augmented sigma points matrix
   */
  MatrixXd AugmentedSigmaPoints();

  /**
   * Create predicted sigma points matrix
   * @param delta_t Time between k and k+1 in s
   */
  void SigmaPointsPrediction(double delta_t, const MatrixXd &Xsig_aug);

  /**
   * Predict radar measurements
   * @param z_out predicted sigma points in measurement space
   * @param S_out predicted measurement covariance matrix
   */
  void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out);

  /**
   * Predict lidar measurements
   * @param z_out predicted measurement
   * @param S_out predicted measurement covariance matrix   
   */
  void PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out);

  /**
   * Common part used by both radar and lidar measurements
   * @param n_z number of measurement states
   * @param Zsig sigma points in measurement space
   * @param z_pred predicted measurement
   */
  MatrixXd PredictMeasurementCommon(const int n_z,
  	                                const MatrixXd &Zsig,
  	                                const VectorXd &z_pred);

  /**
   * Update state mean and covaraince matrix
   * @param n_z number of measurement states
   * @param Zsig sigma points in measurement space
   * @param z_pred predicted measurement
   * @param S predicted measurment covariance matrix
   * @param z ground truth incoming measurement
   */
  void UpdateStateCovariance(const int n_z,
  	                         const MatrixXd &Zsig,
  	                         const VectorXd &z_pred,
	                         const MatrixXd &S,  	                         
  	                         const VectorXd &z);

  /**
   * Consistency check using NIS
   * @param meas_package current sensor measurement
   * @param n_z dimension of sensor measurement
   * @param z_pred prediceted measurement
   * @param S measurement covariance
   * @param e_out output NIS
   */
  void CalculateNIS(MeasurementPackage meas_package,
                    const int n_z,
                    const VectorXd &z_pred,
                    const MatrixXd &S,
                    double* e_out);
};

#endif /* UKF_H */
