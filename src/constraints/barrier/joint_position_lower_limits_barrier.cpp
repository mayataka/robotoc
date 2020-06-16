#include "constraints/barrier/joint_position_lower_limits_barrier.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace barrier {

JointPositionLowerLimits::JointPositionLowerLimits(const Robot& robot, 
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.lowerJointPositionLimit().size()),
    barrier_(barrier),
    qmin_(robot.lowerJointPositionLimit()),
    residual_(Eigen::VectorXd::Zero(qmin_.size())) {
  assert(barrier_ >= 0);
}


void JointPositionLowerLimits::augmentBarrier(const Robot& robot, 
                                              const double dtau, 
                                              const Eigen::VectorXd& q,
                                              Eigen::MatrixXd& Cqq, 
                                              Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  residual_ = dtau * (q-qmin_);
  Cq.array() -= barrier_ / residual_.array();
  for (int i=0; i<dimv_; ++i) {
    Cqq.coeffRef(i, i) += barrier_ / (residual_.coeff(i)*residual_.coeff(i));
  }
}


void JointPositionLowerLimits::augmentBarrierResidual(const Robot& robot, 
                                                      const double dtau, 
                                                      const Eigen::VectorXd& q,
                                                      Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(Cq.size() == robot.dimv());
  residual_ = dtau * (q-qmin_);
  Cq.array() -= barrier_ / residual_.array();
}

} // namespace barrier
} // namespace idocp