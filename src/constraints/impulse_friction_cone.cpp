#include "idocp/constraints/impulse_friction_cone.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

ImpulseFrictionCone::ImpulseFrictionCone(const Robot& robot, 
                                         const double barrier,
                                         const double fraction_to_boundary_rate)
  : ImpulseConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.maxPointContacts()),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
}


ImpulseFrictionCone::ImpulseFrictionCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    fraction_to_boundary_rate_(0) {
}


ImpulseFrictionCone::~ImpulseFrictionCone() {
}


KinematicsLevel ImpulseFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void ImpulseFrictionCone::allocateExtraData(
    ConstraintComponentData& data) const {
  data.r.clear();
  for (int i=0; i<dimc_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(3));
  }
}


bool ImpulseFrictionCone::isFeasible(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      if (frictionConeResidual(mu, s.f[i]) > 0) {
        return false;
      }
    }
  }
  return true;
}


void ImpulseFrictionCone::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    const double mu = robot.frictionCoefficient(i);
    data.slack.coeffRef(i) = - frictionConeResidual(mu, s.f[i]);
  }
  setSlackAndDualPositive(data);
}


void ImpulseFrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      data.r[i].coeffRef(0) = 2 * s.f[i].coeff(0);
      data.r[i].coeffRef(1) = 2 * s.f[i].coeff(1);
      data.r[i].coeffRef(2) = - 2 * mu * mu * s.f[i].coeff(2);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += data.dual.coeff(i) * data.r[i];
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dual_per_slack * data.r[i] * data.r[i].transpose();
      const double coeff 
          = (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
            / data.slack.coeff(i);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias() 
          += coeff * data.r[i];
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s, const ImpulseSplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
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


void ImpulseFrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
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


int ImpulseFrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp