#include "idocp/constraints/contact_normal_force.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

ContactNormalForce::ContactNormalForce(const Robot& robot, const double barrier,
                                       const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.maxPointContacts()) {
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
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      if (s.f[i].coeff(2) < 0) {
        return false;
      }
    }
  }
  return true;
}


void ContactNormalForce::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    data.slack.coeffRef(i) = s.f[i].coeff(2);
  }
  setSlackAndDualPositive(data);
}


void ContactNormalForce::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      kkt_residual.lf().coeffRef(dimf_stack+2) -= dtau * data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void ContactNormalForce::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2) 
          += dtau * data.dual.coeff(i) / data.slack.coeff(i);
      data.residual.coeffRef(i) = - s.f[i].coeff(2) + data.slack.coeff(i);
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
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      data.dslack.coeffRef(i) 
          = d.df().coeff(dimf_stack+2) - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
      dimf_stack += 3;
    }
    else {
      // Set 1.0 to make the fraction-to-boundary rule easy.
      data.slack.coeffRef(i) = 1.0;
      data.dslack.coeffRef(i) = 1.0;
      data.dual.coeffRef(i) = 1.0;
      data.ddual.coeffRef(i) = 1.0;
    }
  }
}


void ContactNormalForce::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      data.residual.coeffRef(i) = - s.f[i].coeff(2) + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
    }
  }
}


int ContactNormalForce::dimc() const {
  return dimc_;
}

} // namespace idocp