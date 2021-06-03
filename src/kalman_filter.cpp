#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_ ;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    // New state update by LASER measurement data
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // convert radar from cartesian to polar coordinates
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    float rho = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    float rho_dot = (px*vx + py*vy) / fmax(rho, 0.001);
    
    VectorXd z_pred = VectorXd(3);
    z_pred << rho, phi, rho_dot;

    VectorXd y = z - z_pred;
    
    while ( y(1) > M_PI || y(1) < -M_PI ) {
        if ( y(1) > M_PI ) {
            y(1) -= 2*M_PI;
        } else {
            y(1) += 2*M_PI;
        }
    }
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    // New state update by Radar measurement data
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
