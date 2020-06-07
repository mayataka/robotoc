#include "ocp/split_ocp.hpp"

#include "Eigen/LU"


namespace invdynocp {

SplitOCP::SplitOCP(const Robot& robot, const CostFunctionInterface* cost, 
                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()) {
  u_.resize(robot.dimv());
  lu_.resize(robot.dimv());
  lx_.resize(2*robot.dimv());
  la_.resize(robot.dimv());
  x_res_.resize(2*robot.dimv());
  lmd_res_.resize(2*robot.dimv());
  a_res_.resize(robot.dimv());
  du_dx_.resize(robot.dimv(), 2*robot.dimv());
  du_da_.resize(robot.dimv(), robot.dimv());
  Qxx_.resize(2*robot.dimv(), 2*robot.dimv());
  Qxa_.resize(2*robot.dimv(), robot.dimv());
  Qaa_.resize(robot.dimv(), robot.dimv());
  Ginv_.resize(robot.dimv(), robot.dimv());
  K_.resize(robot.dimv(), 2*robot.dimv());
  k_.resize(robot.dimv());
}


SplitOCP::~SplitOCP() {
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau,
                            const Eigen::VectorXd& lmd, 
                            const Eigen::VectorXd& x, 
                            const Eigen::VectorXd& a, 
                            const Eigen::VectorXd& lmd_next, 
                            const Eigen::VectorXd& x_next) {
  robot.RNEA(x, a, u_);
  cost_->lu(&robot, t, dtau, x.head(dimq_), x.tail(dimv_), a, u_, lu_);
  cost_->lq(&robot, t, dtau, x.head(dimq_), x.tail(dimv_), a, u_, lx_.head(dimv_));
  cost_->lv(&robot, t, dtau, x.head(dimq_), x.tail(dimv_), a, u_, lx_.tail(dimv_));
  cost_->la(&robot, t, dtau, x.head(dimq_), x.tail(dimv_), a, u_, la_);
  robot.RNEADerivatives(x, a, du_dx_, du_da_);
  lx_.noalias() += du_dx_.transpose() * lu_;
  la_.noalias() += du_da_.transpose() * lu_;
  Qxx_ = lx_ * lx_.transpose();
  Qxa_ = lx_ * la_.transpose();
  Qaa_ = la_ * la_.transpose();
  x_res_ = x - x_next;
  x_res_.head(dimv_).noalias() += dtau * x.tail(dimv_);
  x_res_.tail(dimv_).noalias() += dtau * a;
  lx_.head(dimv_).noalias() += lmd_next.head(dimv_) + lmd_next.tail(dimv_);
  lx_.tail(dimv_).noalias() += dtau * lmd_next.head(dimv_);
  lx_.noalias() -= lmd;
  la_.noalias() += dtau * lmd_next.tail(dimv_);
}


void SplitOCP::linearizeTerminalCost(Robot& robot, const double t, 
                                     const Eigen::VectorXd& x, 
                                     Eigen::MatrixXd& Qxx, 
                                     Eigen::VectorXd& Qx) {
  cost->phiq(&robot, t, x.head(dimq_), x.tail(dimv_), Qx_.head(dimv_));
  cost->phiv(&robot, t, x.head(dimq_), x.tail(dimv_), Qx_.tail(dimv_));
  Qxx = Qx * Qx.transpose();
}


void SplitOCP::backwardRecursion(const double dtau, const Eigen::MatrixXd& P1, 
                                 const Eigen::VectorXd& s1, 
                                 Eigen::MatrixXd& P, Eigen::VectorXd& s) {
  Qxx_.topLeftCorner(dimv_, dimv_).noalias() += P1.topLeftCorner(dimv_, dimv_);
  Qxx_.topLeftCorner(dimv_, dimv_).noalias() += P1.topRightCorner(dimv_, dimv_);
  Qxx_.topLeftCorner(dimv_, dimv_).noalias() += P1.bottomLeftCorner(dimv_, dimv_);
  Qxx_.topLeftCorner(dimv_, dimv_).noalias() += P1.bottomRightCorner(dimv_, dimv_);
  Qxx_.topRightCorner(dimv_, dimv_).noalias() += dtau * P1.topLeftCorner(dimv_, dimv_);
  Qxx_.topRightCorner(dimv_, dimv_).noalias() += dtau * P1.bottomLeftCorner(dimv_, dimv_);
  Qxx_.bottomLeftCorner(dimv_, dimv_).noalias() += dtau * P1.topLeftCorner(dimv_, dimv_);
  Qxx_.bottomLeftCorner(dimv_, dimv_).noalias() += dtau * P1.topRightCorner(dimv_, dimv_);
  const double dtau2 = dtau * dtau;
  Qxx_.bottomRightCorner(dimv_, dimv_).noalias() += dtau2 * P1.topLeftCorner(dimv_, dimv_);
  Qxa_.topRows(dimv_).noalias() += dtau * P1.topRightCorner(dimv_, dimv_);
  Qxa_.topRows(dimv_).noalias() += dtau * P1.bottomRightCorner(dimv_, dimv_);
  Qxa_.bottomRows(dimv_).noalias() += dtau2 * P1.topRightCorner(dimv_, dimv_);
  Qaa_.noalias() += dtau2 * P1.bottomRightCorner(dimv_, dimv_);
  Ginv_ = Qaa_.inverse();
  K_ = - Ginv_ * Qxa_.transpose();
  P = Qxx_;
  P.noalias() += K_.transpose() * Qxa_.transpose();
  lu_ = dtau * s1.tail(dimv_) - la_;
  lu_.noalias() -= dtau * P1.bottomRows(dimv_) * x_res_;
  k_ = Ginv_ * lu_;
  s.head(dimv_) = s1.head(dimv_) + s1.tail(dimv_);
  s.tail(dimv_) = dtau * s1.head(dimv_);
  s.noalias() -= Qxa_ * k_;
  s.noalias() -= lx_;
  lx_ = P.topRows(dimv_) * x_res_;
  s.head(dimv_).noalias() -= lx_;
  s.head(dimv_).noalias() -= P.bottomRows(dimv_) * x_res_;
  s.tail(dimv_).noalias() -= dtau * lx_;
}


void SplitOCP::forwardRecursion(const double dtau, const Eigen::VectorXd& dx, 
                                Eigen::VectorXd& da, 
                                Eigen::VectorXd& dx_next) const {
  du = k_ + K_ * dx;
  dx_next = dx + x_res_;
  dx_next.head(dimv_).noalias() = dtau * dx.tail(dimv_);
  dx_next.tail(dimv_).noalias() = dtau * da;
}


void SplitOCP::updateOCP(Robot& robot, const Eigen::VectorXd& dx, 
                         const Eigen::VectorXd& da, const Eigen::MatrixXd& P, 
                         const Eigen::VectorXd& s, Eigen::VectorXd& x, 
                         Eigen::VectorXd& a, Eigen::VectorXd& lmd) {
  x.noalias() += dx;
  a.noalias() += da;
  lmd.noalias() += P * dx;
  lmd.noalias() -= s;
}


} // namespace invdynocp