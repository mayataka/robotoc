#include "constraints/barrier/joint_velocity_lower_limits_barrier.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace barrier {

JointVelocityLowerLimits::JointVelocityLowerLimits(const Robot& robot,
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    barrier_(barrier),
    vmin_(-robot.jointVelocityLimit()),
    residual_(Eigen::VectorXd::Zero(vmin_.size())) {
  assert(barrier_ >= 0);
  for (int i=0; i<vmin_.size(); ++i) {
    assert(vmin_(i) <= 0);
  }
}


void JointVelocityLowerLimits::augmentBarrier(const Robot& robot, 
                                              const double dtau, 
                                              const Eigen::VectorXd& v,
                                              Eigen::MatrixXd& Cvv, 
                                              Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  residual_ = dtau * (v-vmin_);
  Cv.array() -= barrier_ / residual_.array();
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += barrier_ / (residual_.coeff(i)*residual_.coeff(i));
  }
}


void JointVelocityLowerLimits::augmentBarrierResidual(const Robot& robot, 
                                                      const double dtau, 
                                                      const Eigen::VectorXd& v,
                                                      Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  residual_ = dtau * (v-vmin_);
  Cv.array() -= barrier_ / residual_.array();
}


} // namespace barrier
} // namespace idocp