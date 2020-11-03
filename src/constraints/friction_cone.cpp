#include "idocp/constraints/friction_cone.hpp"

#include <iostream>
#include <assert.h>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, const double barrier,
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(robot.max_point_contacts()) {
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
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      if (frictionConeResidual(mu, s.f[i]) < 0) {
        return false;
      }
    }
  }
  return true;
}


void FrictionCone::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    const double mu = robot.frictionCoefficient(i);
    data.slack.coeffRef(i) = dtau * frictionConeResidual(mu, s.f[i]);
  }
  setSlackAndDualPositive(data);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = 2 * dtau * s.f[i].coeff(0);
      const double fy = 2 * dtau * s.f[i].coeff(1);
      const double fz = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      kkt_residual.lf().coeffRef(dimf_stack  ) += fx * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+1) += fy * data.dual.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack+2) += fz * data.dual.coeff(i);
      dimf_stack += 3;
    }
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = 2 * dtau * s.f[i].coeff(0);
      const double fy = 2 * dtau * s.f[i].coeff(1);
      const double fz = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      const double dual_per_slack = data.dual.coeff(i) / data.slack.coeff(i);
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack  ) 
          += dual_per_slack * fx * fx;
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack+1) 
          += dual_per_slack * fx * fy;
      kkt_matrix.Qff().coeffRef(dimf_stack  , dimf_stack+2) 
          += dual_per_slack * fx * fz;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack  ) 
          += dual_per_slack * fy * fx;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack+1) 
          += dual_per_slack * fy * fy;
      kkt_matrix.Qff().coeffRef(dimf_stack+1, dimf_stack+2) 
          += dual_per_slack * fy * fz;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack  ) 
          += dual_per_slack * fz * fx;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+1) 
          += dual_per_slack * fz * fy;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2) 
          += dual_per_slack * fz * fz;
      data.residual.coeffRef(i) = - dtau * frictionConeResidual(mu, s.f[i]) 
                                  + data.slack.coeff(i);
      data.duality.coeffRef(i) = computeDuality(data.slack.coeff(i), 
                                                data.dual.coeff(i));
      const double coeff 
          = (data.dual.coeff(i)*data.residual.coeff(i)-data.duality.coeff(i)) 
            / data.slack.coeff(i);
      kkt_residual.lf().coeffRef(dimf_stack  ) += fx * coeff;
      kkt_residual.lf().coeffRef(dimf_stack+1) += fy * coeff; 
      kkt_residual.lf().coeffRef(dimf_stack+2) += fz * coeff;
      dimf_stack += 3;
    }
  }
}


void FrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, const SplitDirection& d) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = 2 * dtau * s.f[i].coeff(0);
      const double fy = 2 * dtau * s.f[i].coeff(1);
      const double fz = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      data.dslack.coeffRef(i) = - fx * d.df().coeff(dimf_stack  )
                                - fy * d.df().coeff(dimf_stack+1)
                                - fz * d.df().coeff(dimf_stack+2)
                                - data.residual.coeff(i);
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


void FrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = 2 * dtau * s.f[i].coeff(0);
      const double fy = 2 * dtau * s.f[i].coeff(1);
      const double fz = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      data.residual.coeffRef(i) = - dtau * frictionConeResidual(mu, s.f[i]) 
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