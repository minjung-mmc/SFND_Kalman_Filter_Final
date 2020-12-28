#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  // set weights
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; ++i )
  {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  time_us_ = 0.0;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    // cout << "Kalman Filter Initialization " << endl;

    // set the state with the initial location and zero velocity
    // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad

    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1),
          0.0,
          0.0,
          0.0;
      P_ = MatrixXd::Identity(n_x_, n_x_);
      // P_(0,0) = 0.5; P_(1, 1) = 0.5;
      P_(2,2) = 1.5; P_(3,3) = 1.5; P_(4, 4) = 1.5;

      // std::cout << "<< " << x_ << ">>" << std::endl;
      std::cout << "Initialized by Lidar" << std::endl;
    }
    else
    {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      double px = rho * cos(phi);
      double py = rho * sin(phi);

      x_ << px,
          py,
          0.0,
          0.0,
          0.0;
      P_ = MatrixXd::Identity(n_x_, n_x_);
      std::cout << "Initialized by Radar" << std::endl;
    }
    is_initialized_ = true;
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; // microseconds to seconds

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }

  time_us_ = meas_package.timestamp_;
}


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out )
{
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug << x_, 0, 0; 
  Xsig_aug.col(0) = x_aug;
  // create augmented covariance matrix
  MatrixXd nu = MatrixXd(2, 2);
  nu << std_a_, 0, 0, std_yawdd_;
  MatrixXd Q = nu.transpose() * nu;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  // write result
  *Xsig_out = Xsig_aug;  
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t) {
    // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  MatrixXd Xsig_aug = *Xsig_out;

  // predict sigma points
  VectorXd deltaX = VectorXd(n_x_);
  VectorXd noise = VectorXd(n_x_);
  VectorXd X = VectorXd(n_x_);
  double v, yaw, yawRate, std_a, std_yawdd;

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    v = Xsig_aug(2, i);
    yaw = Xsig_aug(3, i);
    yawRate = Xsig_aug(4, i);
    std_a = Xsig_aug(5, i);
    std_yawdd = Xsig_aug(6, i);

    if (fabs(yawRate) > 1e-3)
    {
      deltaX << v / yawRate * (sin(yaw + yawRate * delta_t) - sin(yaw)), v / yawRate * (-cos(yaw + yawRate * delta_t) + cos(yaw)), 0, yawRate * delta_t, 0;
    }
    else
    {
      deltaX << v * cos(yaw) * delta_t, v * sin(yaw) * delta_t, 0, 0, 0;
    }
    noise << 0.5 * pow(delta_t, 2) * cos(yaw) * std_a, 
    0.5 * pow(delta_t, 2) * sin(yaw) * std_a, 
    delta_t * std_a, 
    0.5 * pow(delta_t, 2) * std_yawdd, 
    delta_t * std_yawdd;

    X << Xsig_aug(0, i), Xsig_aug(1, i), v, yaw, yawRate;
    Xsig_pred.col(i) = X + deltaX + noise;
    
  }
  // avoid division by zero
  *Xsig_out = Xsig_pred;
}


void UKF::PredictMeanAndCovariance(MatrixXd* Xsig_out) {

MatrixXd Xsig_pred = *Xsig_out;
VectorXd x = VectorXd(n_x_);
MatrixXd P = MatrixXd(n_x_, n_x_);                  
// predict state mean
for (int i = 0; i < 2 * n_aug_ + 1; ++i)
{
  x += weights_(i) * Xsig_pred.col(i);
}

// predict state covariance matrix
P.fill(0.0);
for (int i = 0; i < 2 * n_aug_ + 1; ++i)
{
  VectorXd x_diff = Xsig_pred.col(i) - x;
  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
  P += weights_(i) * (x_diff) * (x_diff).transpose();
}

x_ = x;
P_ = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Xsig_out, MatrixXd* Zsig_out) {
MatrixXd Xsig_pred = *Xsig_out;
int n_z_ = 3;

// create matrix for sigma points in measurement space
MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

// mean predicted measurement
VectorXd z_pred = VectorXd(n_z_);
z_pred.fill(0.0);
// measurement covariance matrix S
MatrixXd S = MatrixXd(n_z_, n_z_);
S.fill(0.0);

 // transform sigma points into measurement space
  double px, py, v, yaw, yawRate, rho, phi, rhoDot;
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      px = Xsig_pred(0,i); py = Xsig_pred(1,i); v = Xsig_pred(2,i); yaw = Xsig_pred(3,i); yawRate = Xsig_pred(4,i);
      rho = sqrt(pow(px,2)+pow(py,2)); 
      phi = atan2(py,px);
      rhoDot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
      Zsig.col(i) << rho, phi, rhoDot;
      z_pred += weights_(i)*Zsig.col(i); 
  }

  for (int i = 0; i< 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S += weights_(i)*(z_diff)*(z_diff).transpose();
  }
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_,0,0,0,std_radphi_,0,0,0,std_radrd_;
  R = R*R;
  S += R;
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

// Used in UKF results
void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // x_, P_, Xsig_pred_ is changed.
  MatrixXd Xsig_out = MatrixXd(n_aug_, 2*n_aug_+1);
  AugmentedSigmaPoints(&Xsig_out);
  SigmaPointPrediction(&Xsig_out, delta_t);
  PredictMeanAndCovariance(&Xsig_out);
  Xsig_pred_ = Xsig_out;
}


void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z_ = 2;

    VectorXd z(n_z_);
    z << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1];

    // measurement covariance
    MatrixXd R_ = MatrixXd(2, 2);
    R_ << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;

    // measurement matrix
    MatrixXd H_ = MatrixXd(2, 5);
    H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    //new estimate
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  int n_z_ = 3;

  VectorXd z(n_z_);
  z << meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1],
      meas_package.raw_measurements_[2];

  MatrixXd Zsig_out = MatrixXd(n_z_, 2*n_aug_+1);
  VectorXd z_out = VectorXd(n_z_);
  MatrixXd S_out = MatrixXd(n_z_, n_z_);
  PredictRadarMeasurement(&z_out, &S_out, &Xsig_pred_, &Zsig_out);
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  // calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_ +1; ++i)
  {
    VectorXd z_diff = Zsig_out.col(i) - z_out;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      Tc += weights_(i)*(x_diff)*(z_diff).transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z_);
  K = Tc * S_out.inverse();

  VectorXd z_diff = z - z_out;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S_out*K.transpose();
}