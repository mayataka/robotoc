#include "constraints/joint_space_constraints.hpp"

#include <assert.h>


namespace idocp {


JointSpaceConstraints::JointSpaceConstraints(const Robot& robot) 
  : q_max_(robot.upperJointPositionLimit()),
    q_min_(robot.lowerJointPositionLimit()),
    v_max_(robot.jointVelocityLimit()),
    v_min_(-robot.jointVelocityLimit()),
    u_max_(robot.jointEffortLimit()),
    u_min_(-robot.jointEffortLimit()),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()) {
}


void JointSpaceConstraints::Cq(const Robot& robot, const double dtau, 
                               const Eigen::VectorXd& q, Eigen::VectorXd& Cq) {
  assert(q.size() == q_max.size());
  assert(Cq.size() == 2*q_max.size());
  Cq.head(dimq_) = tau * (q_max-q);
  Cq.tail(dimq_) = tau * (q-q_min);
}


void JointSpaceConstraints::Cv(const Robot& robot, const double dtau, 
                               const Eigen::VectorXd& v, Eigen::VectorXd& Cv) {
  assert(v.size() == v_max.size());
  assert(Cv.size() == 2*v_max.size());
  Cv.head(dimv_) = tau * (v_max-v);
  Cv.tail(dimv_) = tau * (v-v_min);
}


void JointSpaceConstraints::Cu(const Robot& robot, const double dtau, 
                               const Eigen::VectorXd& u, Eigen::VectorXd& Cu) {
  assert(u.size() == u_max.size());
  assert(Cu.size() == 2*u_max.size());
  Cv.head(dimv_) = tau * (u_max-u);
  Cv.tail(dimv_) = tau * (u-u_min);
}

} // namespace idocp