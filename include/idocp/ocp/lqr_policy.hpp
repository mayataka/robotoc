#ifndef IDOCP_LQR_POLICY_HPP_
#define IDOCP_LQR_POLICY_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

struct LQRPolicy {
  LQRPolicy(const Robot& robot)
    : K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
      k(Eigen::VectorXd::Zero(robot.dimu())) {
  }

  LQRPolicy() 
    : K(),
      k() {
  }

  ~LQRPolicy() {
  }

  Eigen::MatrixXd K;
  Eigen::VectorXd k;

};

} // namespace idocp 

#endif // IDOCP_LQR_POLICY_HPP_ 