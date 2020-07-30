#include "mpc_simulator/runge_kutta.hpp"

#include <assert.h>


namespace idocp {
namespace simulator {

RungeKutta::RungeKutta(const Robot& robot) 
  : robot_(robot),
    dq_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_(Eigen::VectorXd::Zero(robot.dimv())),
    kq1_(Eigen::VectorXd::Zero(robot.dimq())),
    kq2_(Eigen::VectorXd::Zero(robot.dimq())),
    kq3_(Eigen::VectorXd::Zero(robot.dimq())),
    kq4_(Eigen::VectorXd::Zero(robot.dimq())),
    kv1_(Eigen::VectorXd::Zero(robot.dimv())),
    kv2_(Eigen::VectorXd::Zero(robot.dimv())),
    kv3_(Eigen::VectorXd::Zero(robot.dimv())),
    kv4_(Eigen::VectorXd::Zero(robot.dimv())) {
}


void RungeKutta::integrate(const double integration_length, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& tau, Eigen::VectorXd& q_next, 
                           Eigen::VectorXd& v_next) {
  assert(integration_length > 0);
  assert(q.size() == robot_.dimq());
  assert(v.size() == robot_.dimv());
  assert(tau.size() == robot_.dimv());
  assert(q_next.size() == robot_.dimq());
  assert(v_next.size() == robot_.dimv());
  robot_.stateEquation(q, v, tau, dq_, dv_);
  kq1_ = integration_length * dq_;
  kv1_ = integration_length * dv_;
  robot_.stateEquation(q+0.5*kq1_, v+0.5*kv1_, tau, dq_, dv_);
  kq2_ = integration_length * dq_;
  kv2_ = integration_length * dv_;
  robot_.stateEquation(q+0.5*kq2_, v+0.5*kv2_, tau, dq_, dv_);
  kq3_ = integration_length * dq_;
  kv3_ = integration_length * dv_;
  robot_.stateEquation(q+kq3_, v+kv3_, tau, dq_, dv_);
  kq4_ = integration_length * dq_;
  kv4_ = integration_length * dv_;
  q_next = q + (kq1_+2*kq2_+2*kq3_+kq4_) / 6;
  v_next = v + (kv1_+2*kv2_+2*kv3_+kv4_) / 6;
}

} // namespace simulator
} // namespace idocp
