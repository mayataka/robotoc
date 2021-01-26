#include "idocp/constraints/friction_cone.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, const double barrier,
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.maxPointContacts()),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
}


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    fraction_to_boundary_rate_(0) {
}


FrictionCone::~FrictionCone() {
}


bool FrictionCone::useKinematics() const {
  return false;
}


KinematicsLevel FrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void FrictionCone::allocateExtraData(ConstraintComponentData& data) const {
  data.r.clear();
  for (int i=0; i<dimc_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(3));
  }
}


bool FrictionCone::isFeasible(Robot& robot, ConstraintComponentData& data, 
                              const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      if (frictionConeResidual(mu, s.f[i]) > 0) {
        return false;
      }
    }
  }
  return true;
}


void FrictionCone::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    const double mu = robot.frictionCoefficient(i);
    data.slack.coeffRef(i) = - frictionConeResidual(mu, s.f[i]);
  }
  setSlackAndDualPositive(data);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      data.r[i].coeffRef(0) = 2 * s.f[i].coeff(0);
      data.r[i].coeffRef(1) = 2 * s.f[i].coeff(1);
      data.r[i].coeffRef(2) = - 2 * mu * mu * s.f[i].coeff(2);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dtau * data.dual.coeff(i) * data.r[i];
      dimf_stack += 3;
    }
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dtau * dual_per_slack * data.r[i] * data.r[i].transpose();
      const double coeff 
          = (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
            / data.slack.coeff(i);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias() 
          += dtau * coeff * data.r[i];
      dimf_stack += 3;
    }
  }
}


void FrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      data.dslack.coeffRef(i) = - data.r[i].dot(d.df().template segment<3>(dimf_stack))
                                - data.residual.coeff(i);
      data.ddual.coeffRef(i) = computeDualDirection(data.slack.coeff(i), 
                                                    data.dual.coeff(i), 
                                                    data.dslack.coeff(i), 
                                                    data.duality.coeff(i));
      dimf_stack += 3;
    }
    else {
      data.residual.coeffRef(i) = 0;
      data.duality.coeffRef(i)  = 0;
      data.slack.coeffRef(i)    = 1.0;
      data.dslack.coeffRef(i)   = fraction_to_boundary_rate_;
      data.dual.coeffRef(i)     = 1.0;
      data.ddual.coeffRef(i)    = fraction_to_boundary_rate_;
    }
  }
}


void FrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      data.residual.coeffRef(i) = frictionConeResidual(mu, s.f[i]) 
                                  + data.slack.coeff(i);
      data.duality.coeffRef(i)  = computeDuality(data.slack.coeff(i), 
                                                 data.dual.coeff(i));
    }
    else {
      data.residual.coeffRef(i) = 0;
      data.duality.coeffRef(i)  = 0;
    }
  }
}


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp