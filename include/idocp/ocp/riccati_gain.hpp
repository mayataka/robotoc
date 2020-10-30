#ifndef IDOCP_RICCATI_GAIN_HPP_
#define IDOCP_RICCATI_GAIN_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiGain {
public:
  RiccatiGain(const Robot& robot)
    : K(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())),
      k(Eigen::VectorXd::Zero(robot.dimu())),
      dimv_(robot.dimv()),
      dimu_(robot.dimu()) {
  }

  RiccatiGain()
    : K(),
      k(),
      dimv_(0),
      dimu_(0) {
  }

  ~RiccatiGain() {
  }

  RiccatiGain(const RiccatiGain&) = default;

  RiccatiGain& operator=(const RiccatiGain&) = default;
 
  RiccatiGain(RiccatiGain&&) noexcept = default;

  RiccatiGain& operator=(RiccatiGain&&) noexcept = default;

  const Eigen::Block<const Eigen::MatrixXd> Kq() const {
    return K.topLeftCorner(dimu_, dimv_);
  }

  const Eigen::Block<const Eigen::MatrixXd> Kv() const {
    return K.topRightCorner(dimu_, dimv_);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::MatrixXd K;
  Eigen::VectorXd k;

private:
  int dimv_, dimu_;

};

} // namespace idocp 

#endif // IDOCP_RICCATI_GAIN_HPP_