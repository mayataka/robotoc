#include "idocp/constraints/contact_normal_force.hpp"

#include <iostream>
#include <assert.h>


namespace idocp {

ContactNormalForce::ContactNormalForce(const Robot& robot, const double barrier,
                                       const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.max_point_contacts()) {
}


ContactNormalForce::ContactNormalForce()
  : ConstraintComponentBase(),
    dimc_(0) {
}


ContactNormalForce::~ContactNormalForce() {
}


bool ContactNormalForce::useKinematics() const {
  return false;
}


KinematicsLevel ContactNormalForce::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


bool ContactNormalForce::isFeasible(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s) const {
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      if (s.f[i].coeff(2) < 0) {
        return false;
      }
    }
  }
  return true;
}


void ContactNormalForce::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    const double mu = robot.frictionCoefficient(i);
    data.slack.coeffRef(i) = dtau * s.f[i].coeff(2);
  }
  setSlackAndDualPositive(data);
}


void ContactNormalForce::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      kkt_residual.lf().coeffRef(dimf_stack+2) -= dtau * data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void ContactNormalForce::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
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


void ContactNormalForce::computeSlackAndDualDirection(
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


void ContactNormalForce::computePrimalAndDualResidual(
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


int ContactNormalForce::dimc() const {
  return dimc_;
}

} // namespace idocp