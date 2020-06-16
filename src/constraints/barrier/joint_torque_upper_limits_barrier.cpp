#include "constraints/barrier/joint_torque_upper_limits_barrier.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace barrier {

JointTorqueUpperLimits::JointTorqueUpperLimits(const Robot& robot,
                                               const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointEffortLimit().size()),
    barrier_(barrier),
    umax_(robot.jointEffortLimit()),
    residual_(Eigen::VectorXd::Zero(umax_.size())) {
  assert(barrier_ >= 0);
  for (int i=0; i<umax_.size(); ++i) {
    assert(umax_(i) >= 0);
  }
}


void JointTorqueUpperLimits::augmentBarrier(const Robot& robot, 
                                            const double dtau, 
                                            const Eigen::VectorXd& u,
                                            Eigen::MatrixXd& Cuu, 
                                            Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(Cuu.rows() == robot.dimv());
  assert(Cuu.cols() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  residual_ = dtau * (umax_-u);
  Cu.array() -= barrier_ / residual_.array();
  for (int i=0; i<dimv_; ++i) {
    Cuu.coeffRef(i, i) += barrier_ / (residual_.coeff(i)*residual_.coeff(i));
  }
}


void JointTorqueUpperLimits::augmentBarrierResidual(const Robot& robot, 
                                                    const double dtau, 
                                                    const Eigen::VectorXd& u,
                                                    Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  residual_ = dtau * (umax_-u);
  Cu.array() -= barrier_ / residual_.array();
}


} // namespace barrier
} // namespace idocp