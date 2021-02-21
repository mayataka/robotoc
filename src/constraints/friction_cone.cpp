#include "idocp/constraints/friction_cone.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

FrictionCone::FrictionCone(const Robot& robot, const double mu, 
                           const double barrier,
                           const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
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


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    mu_(0) {
}


FrictionCone::~FrictionCone() {
}


void FrictionCone::setFrictionCoefficient(const double mu) {
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


void FrictionCone::setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                                   const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    const int idx1 = 2*i;
    data.slack.coeffRef(idx1) = - normalForceResidual(s.f[i]);
    const int idx2 = 2*i+1;
    data.slack.coeffRef(idx2) = - frictionConeResidual(mu_, s.f[i]);
  }
  setSlackAndDualPositive(data);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx1 = 2*i;
      kkt_residual.lf().coeffRef(dimf_stack+2) -= dt * data.dual.coeff(idx1);
      const int idx2 = 2*i+1;
      data.r[i].coeffRef(0) = 2 * s.f[i].coeff(0);
      data.r[i].coeffRef(1) = 2 * s.f[i].coeff(1);
      data.r[i].coeffRef(2) = - 2 * mu_ * mu_ * s.f[i].coeff(2);
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dt * data.dual.coeff(idx2) * data.r[i];
      dimf_stack += 3;
    }
  }
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx1 = 2*i;
      kkt_matrix.Qff().coeffRef(dimf_stack+2, dimf_stack+2)
          += dt * (data.dual.coeff(idx1) / data.slack.coeff(idx1));
      kkt_residual.lf().coeffRef(dimf_stack+2) 
          -= dt * (data.dual.coeff(idx1)*data.residual.coeff(idx1)-data.duality.coeff(idx1)) 
              / data.slack.coeff(idx1);
      const int idx2 = 2*i+1;
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dt * (data.dual.coeff(idx2) / data.slack.coeff(idx2)) 
                  * data.r[i] * data.r[i].transpose();
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dt * (data.dual.coeff(idx2)*data.residual.coeff(idx2)-data.duality.coeff(idx2)) 
              / data.slack.coeff(idx2)  * data.r[i];
      dimf_stack += 3;
    }
  }
}


void FrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  // Because data.slack(i) and data.dual(i) are always positive,  
  // - fraction_rate * (slack.coeff(i)/dslack.coeff(i)) and 
  // - fraction_rate * (dual.coeff(i)/ddual.coeff(i))  
  // at the inactive constraint index i are always negative, 
  // and therefore do not affect to step size.
  data.dslack.fill(1.0);
  data.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
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


void FrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual.setZero();
  data.duality.setZero();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
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


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp