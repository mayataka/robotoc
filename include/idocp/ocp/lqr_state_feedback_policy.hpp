#ifndef IDOCP_LQR_STATE_FEEDBACK_POLICY_HPP_
#define IDOCP_LQR_STATE_FEEDBACK_POLICY_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class LQRStateFeedbackPolicy {
public:
  LQRStateFeedbackPolicy(const Robot& robot)
    : K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
      k(Eigen::VectorXd::Zero(robot.dimu())) {
  }

  LQRStateFeedbackPolicy() 
    : K(),
      k() {
  }

  ~LQRStateFeedbackPolicy() {
  }

  LQRStateFeedbackPolicy(const LQRStateFeedbackPolicy&) = default;

  LQRStateFeedbackPolicy& operator=(const LQRStateFeedbackPolicy&) = default;
 
  LQRStateFeedbackPolicy(LQRStateFeedbackPolicy&&) noexcept = default;

  LQRStateFeedbackPolicy& operator=(LQRStateFeedbackPolicy&&) noexcept = default;

  Eigen::MatrixXd K;
  Eigen::VectorXd k;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

} // namespace idocp 

#endif // IDOCP_LQR_STATE_FEEDBACK_POLICY_HPP_ 