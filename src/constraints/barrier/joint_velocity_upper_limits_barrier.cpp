#include "constraints/barrier/joint_velocity_upper_limits_barrier.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace barrier {

JointVelocityUpperLimits::JointVelocityUpperLimits(const Robot& robot,
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    barrier_(barrier),
    vmax_(robot.jointVelocityLimit()),
    residual_(Eigen::VectorXd::Zero(vmax_.size())) {
  assert(barrier_ >= 0);
  for (int i=0; i<vmax_.size(); ++i) {
    assert(vmax_(i) >= 0);
  }
}


void JointVelocityUpperLimits::augmentBarrier(const Robot& robot, 
                                              const double dtau, 
                                              const Eigen::VectorXd& v,
                                              Eigen::MatrixXd& Cvv, 
                                              Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  residual_ = dtau * (vmax_-v);
  Cv.array() -= barrier_ / residual_.array();
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += barrier_ / (residual_.coeff(i)*residual_.coeff(i));
  }
}


void JointVelocityUpperLimits::augmentBarrierResidual(const Robot& robot, 
                                                      const double dtau, 
                                                      const Eigen::VectorXd& v,
                                                      Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  residual_ = dtau * (vmax_-v);
  Cv.array() -= barrier_ / residual_.array();
}


} // namespace barrier
} // namespace idocp