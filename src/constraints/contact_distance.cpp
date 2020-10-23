#include "idocp/constraints/contact_distance.hpp"

#include <iostream>
#include <assert.h>


namespace idocp {

ContactDistance::ContactDistance(const Robot& robot, const double barrier,
                                 const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.max_point_contacts()) {
}


ContactDistance::ContactDistance()
  : ConstraintComponentBase(),
    dimc_(0) {
}


ContactDistance::~ContactDistance() {
}


bool ContactDistance::useKinematics() const {
  return false;
}


KinematicsLevel ContactDistance::kinematicsLevel() const {
  return KinematicsLevel::PositionLevel;
}


bool ContactDistance::isFeasible(Robot& robot, ConstraintComponentData& data, 
                                 const SplitSolution& s) const {
  Eigen::Vector3d end_effector_position = Eigen::Vector3d::Zero();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!s.isContactActive(i)) {
      robot.computeContactResidual(i, end_effector_position);
      if (end_effector_position.coeff(2) < 0) {
        return false;
      }
    }
  }
  return true;
}


void ContactDistance::setSlackAndDual(Robot& robot, 
                                      ConstraintComponentData& data, 
                                      const double dtau, 
                                      const SplitSolution& s) const {
  assert(dtau > 0);
  Eigen::Vector3d end_effector_position = Eigen::Vector3d::Zero();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    robot.computeContactResidual(i, end_effector_position);
    data.slack.coeffRef(i) = dtau * end_effector_position.coeff(2);
  }
  setSlackAndDualPositive(data);
}


void ContactDistance::augmentDualResidual(Robot& robot, 
                                          ConstraintComponentData& data, 
                                          const double dtau, 
                                          const SplitSolution& s, 
                                          KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  Eigen::MatrixXd end_effector_Jacobian = Eigen::MatrixXd::Zero(s.dimf(), 
                                                                robot.dimv());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!s.isContactActive(i)) {
      robot.computeContactDerivative(i, end_effector_Jacobian);
    }
  }
  if (s.dimf() > 0) {
    kkt_residual.lq().noalias -= dtau * end_effector_Jacobian.row(2) 
                                      * data.dual.coeff(i);
  }
}


void ContactDistance::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  Eigen::MatrixXd end_effector_Jacobian = Eigen::MatrixXd::Zero(3, robot.dimv());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (!robot.is_contact_active(i)) {
      robot.computeContactDerivative(i, end_effector_Jacobian);
      kkt_matrix.Qqq().noalias() 
          += dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i) 
                         * frame_derivative_[i].row(2).transpose() 
                         * frame_derivative_[i].row(2);
      robot.computeContactResidual(i, frame_position_[i]);
      data.residual.coeffRef(i) = - dtau * frame_position_[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      kkt_residual.lq().noalias()
          -= dtau * frame_derivative_[i].row(2)
              * (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
                  / data.slack.coeff(i);
    }
  }




  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2) 
          = dtau * dtau * data.dual.coeff(i) / data.slack.coeff(i);
      data.residual.coeffRef(i) = - dtau * s.f[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      kkt_residual.lf().coeffRef(dimf_stack+2) 
          -= dtau * (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
                  / data.slack.coeff(i);
      dimf_stack += 3;
    }
  }
}


void ContactDistance::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, const SplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      data.dslack.coeffRef(i) 
          = dtau * d.df().coeff(dimf_stack+2) - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
      dimf_stack += 3;
    }
    else {
      data.slack.coeffRef(i) = 1;
      data.dslack.coeffRef(i) = 1;
      data.dual.coeffRef(i) = 1;
      data.ddual.coeffRef(i) = 1;
    }
  }
}


void ContactDistance::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      data.residual.coeffRef(i) = - dtau * s.f[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
    }
  }
}


int ContactDistance::dimc() const {
  return dimc_;
}

} // namespace idocp