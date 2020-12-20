#include "idocp/constraints/friction_cone.hpp"

#include <iostream>
#include <cassert>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, const double barrier,
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.maxPointContacts()) {
}


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0) {
}


FrictionCone::~FrictionCone() {
}


bool FrictionCone::useKinematics() const {
  return false;
}


KinematicsLevel FrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
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
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double gx = 2 * dtau * s.f[i].coeff(0);
      const double gy = 2 * dtau * s.f[i].coeff(1);
      const double gz = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      kkt_residual.lf().coeffRef(dimf_stack  ) += gx * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+1) += gy * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+2) += gz * data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double gx = 2 * s.f[i].coeff(0);
      const double gy = 2 * s.f[i].coeff(1);
      const double gz = - 2 * mu * mu * s.f[i].coeff(2);
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      const double dtau_dual_per_slack = dtau * dual_per_slack;
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack  ) += dtau_dual_per_slack * gx * gx;
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack+1) += dtau_dual_per_slack * gx * gy;
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack+2) += dtau_dual_per_slack * gx * gz;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack  ) += dtau_dual_per_slack * gy * gx;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack+1) += dtau_dual_per_slack * gy * gy;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack+2) += dtau_dual_per_slack * gy * gz;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack  ) += dtau_dual_per_slack * gz * gx;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+1) += dtau_dual_per_slack * gz * gy;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2) += dtau_dual_per_slack * gz * gz;
      data.residual.coeffRef(i) = frictionConeResidual(mu, s.f[i]) 
                                  + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      const double coeff 
          = (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
            / data.slack.coeff(i);
      const double dtau_coeff = dtau * coeff;
      kkt_residual.lf().coeffRef(dimf_stack  ) += gx * dtau_coeff;
      kkt_residual.lf().coeffRef(dimf_stack+1) += gy * dtau_coeff; 
      kkt_residual.lf().coeffRef(dimf_stack+2) += gz * dtau_coeff;
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
      const double mu = robot.frictionCoefficient(i);
      const double gx = 2 * s.f[i].coeff(0);
      const double gy = 2 * s.f[i].coeff(1);
      const double gz = - 2 * mu * mu * s.f[i].coeff(2);
      data.dslack.coeffRef(i) = - gx * d.df().coeff(dimf_stack  )
                                - gy * d.df().coeff(dimf_stack+1)
                                - gz * d.df().coeff(dimf_stack+2)
                                - data.residual.coeff(i);
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


void FrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      data.residual.coeffRef(i) = frictionConeResidual(mu, s.f[i]) 
                                  + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
    }
  }
}


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp