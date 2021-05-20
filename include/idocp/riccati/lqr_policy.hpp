#ifndef IDOCP_LQR_POLICY_HPP_
#define IDOCP_LQR_POLICY_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class LQRPolicy {
public:
  LQRPolicy(const Robot& robot)
    : K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
      k(Eigen::VectorXd::Zero(robot.dimu())),
      dimv_(robot.dimv()),
      dimu_(robot.dimu()) {
  }

  LQRPolicy() 
    : K(),
      k(),
      dimv_(0),
      dimu_(0) {
  }

  ~LQRPolicy() {
  }

  Eigen::MatrixXd K;

  Eigen::VectorXd k;

  const Eigen::Block<const Eigen::MatrixXd> Kq() const {
    return K.topLeftCorner(dimu_, dimv_);
  }

  const Eigen::Block<const Eigen::MatrixXd> Kv() const {
    return K.topRightCorner(dimu_, dimv_);
  }

  bool isApprox(const LQRPolicy& other) const {
    if (!K.isApprox(other.K)) return false;
    if (!k.isApprox(other.k)) return false;
    return true;
  }

private:
  int dimv_, dimu_;

};

} // namespace idocp 

#endif // IDOCP_LQR_POLICY_HPP_ 