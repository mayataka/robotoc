#include "idocp/constraints/impulse_friction_cone.hpp"

#include <iostream>


namespace idocp {

ImpulseFrictionCone::ImpulseFrictionCone(const Robot& robot, const double mu,
                                         const double barrier,
                                         const double fraction_to_boundary_rate)
  : ImpulseConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(2*robot.maxPointContacts()),
    mu_(mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range("invalid value: mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


ImpulseFrictionCone::ImpulseFrictionCone()
  : ImpulseConstraintComponentBase(),
    dimc_(0),
    mu_(0) {
}


ImpulseFrictionCone::~ImpulseFrictionCone() {
}


void ImpulseFrictionCone::setFrictionCoefficient(const double mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range("invalid value: mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  mu_ = mu;
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


bool ImpulseFrictionCone::isFeasible(Robot& robot, 
                                     ConstraintComponentData& data, 
                                     const ImpulseSplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      if (frictionConeResidual(mu_, s.f[i]) > 0) {
        return false;
      }
      if (normalForceResidual(s.f[i]) > 0) {
        return false;
      }
    }
  }
  return true;
}


void ImpulseFrictionCone::setSlackAndDual(Robot& robot, 
                                          ConstraintComponentData& data, 
                                          const ImpulseSplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    const int idx1 = 2*i;
    data.slack.coeffRef(idx1) = - normalForceResidual(s.f[i]);
    const int idx2 = 2*i+1;
    data.slack.coeffRef(idx2) = - frictionConeResidual(mu_, s.f[i]);
  }
  setSlackAndDualPositive(data);
}


void ImpulseFrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx1 = 2*i;
      kkt_residual.lf().coeffRef(dimf_stack+2) -= data.dual.coeff(idx1);
      const int idx2 = 2*i+1;
      data.r[i].coeffRef(0) = 2 * s.f[i].coeff(0);
      data.r[i].coeffRef(1) = 2 * s.f[i].coeff(1);
      data.r[i].coeffRef(2) = - 2 * mu_ * mu_ * s.f[i].coeff(2);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += data.dual.coeff(idx2) * data.r[i];
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx1 = 2*i;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2)
          += (data.dual.coeff(idx1) / data.slack.coeff(idx1));
      kkt_residual.lf().coeffRef(dimf_stack+2) 
          -= (data.dual.coeff(idx1)*data.residual.coeff(idx1)-data.duality.coeff(idx1)) 
              / data.slack.coeff(idx1);
      const int idx2 = 2*i+1;
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += (data.dual.coeff(idx2) / data.slack.coeff(idx2)) 
              * data.r[i] * data.r[i].transpose();
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += (data.dual.coeff(idx2)*data.residual.coeff(idx2)-data.duality.coeff(idx2)) 
              / data.slack.coeff(idx2)  * data.r[i];
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const ImpulseSplitSolution& s, 
    const ImpulseSplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // - fraction_rate * (slack.coeff(i)/dslack.coeff(i)) and 
  // - fraction_rate * (dual.coeff(i)/ddual.coeff(i))  
  // at the inactive constraint index i are always negative, 
  // and therefore do not affect to step size.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx1 = 2*i;
      data.dslack.coeffRef(idx1) = d.df().coeff(dimf_stack+2) 
                                    - data.residual.coeff(idx1);
      data.ddual.coeffRef(idx1)  = computeDualDirection(data.slack.coeff(idx1), 
                                                        data.dual.coeff(idx1), 
                                                        data.dslack.coeff(idx1), 
                                                        data.duality.coeff(idx1));
      const int idx2 = 2*i+1;
      data.dslack.coeffRef(idx2) 
          = - data.r[i].dot(d.df().template segment<3>(dimf_stack)) 
            - data.residual.coeff(idx2);
      data.ddual.coeffRef(idx2) = computeDualDirection(data.slack.coeff(idx2), 
                                                       data.dual.coeff(idx2), 
                                                       data.dslack.coeff(idx2), 
                                                       data.duality.coeff(idx2));
      dimf_stack += 3;
    }
  }
}


void ImpulseFrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, 
    const ImpulseSplitSolution& s) const {
  data.residual.setZero();
  data.duality.setZero();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      const int idx1 = 2*i;
      data.residual.coeffRef(idx1) = normalForceResidual(s.f[i]) 
                                      + data.slack.coeff(idx1);
      data.duality.coeffRef(idx1)  = computeDuality(data.slack.coeff(idx1), 
                                                    data.dual.coeff(idx1));
      const int idx2 = 2*i+1;
      data.residual.coeffRef(idx2) = frictionConeResidual(mu_, s.f[i]) 
                                      + data.slack.coeff(idx2);
      data.duality.coeffRef(idx2)  = computeDuality(data.slack.coeff(idx2), 
                                                    data.dual.coeff(idx2));
    }
  }
}


int ImpulseFrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp