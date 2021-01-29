#include "idocp/constraints/linearized_friction_cone.hpp"

#include <iostream>
#include <cassert>
#include <stdexcept>


namespace idocp {

LinearizedFrictionCone::LinearizedFrictionCone(
    const Robot& robot, const double mu, const double barrier,
    const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    dimc_(5*robot.maxPointContacts()),
    mu_(mu),
    Jac_() {
  try {
    if (mu <= 0) {
      throw std::out_of_range("invalid value: mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  Jac_ <<  0,  0, -1, 
           1,  0, -(mu/std::sqrt(2)),
          -1,  0, -(mu/std::sqrt(2)),
           0,  1, -(mu/std::sqrt(2)),
           0, -1, -(mu/std::sqrt(2));
}


LinearizedFrictionCone::LinearizedFrictionCone()
  : ConstraintComponentBase(),
    dimc_(0),
    mu_(0),
    Jac_() {
}


LinearizedFrictionCone::~LinearizedFrictionCone() {
}


void LinearizedFrictionCone::setFrictionCoefficient(const double mu) {
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
  Jac_ <<  0,  0, -1, 
           1,  0, -(mu/std::sqrt(2)),
          -1,  0, -(mu/std::sqrt(2)),
           0,  1, -(mu/std::sqrt(2)),
           0, -1, -(mu/std::sqrt(2));
}


bool LinearizedFrictionCone::useKinematics() const {
  return false;
}


KinematicsLevel LinearizedFrictionCone::kinematicsLevel() const {
  return KinematicsLevel::AccelerationLevel;
}


void LinearizedFrictionCone::allocateExtraData(
    ConstraintComponentData& data) const {
  data.r.clear();
  for (int i=0; i<dimc_; ++i) {
    data.r.push_back(Eigen::VectorXd::Zero(5));
  }
  data.J.clear();
  for (int i=0; i<dimc_; ++i) {
    data.J.push_back(Eigen::MatrixXd::Zero(5, 3));
  }
}


bool LinearizedFrictionCone::isFeasible(Robot& robot, 
                                        ConstraintComponentData& data, 
                                        const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      frictionConeResidual(mu_, s.f[i], data.residual.template segment<5>(5*i));
      if (data.residual.maxCoeff() > 0) {
        return false;
      }
    }
  }
  return true;
}


void LinearizedFrictionCone::setSlackAndDual(Robot& robot, 
                                             ConstraintComponentData& data, 
                                             const SplitSolution& s) const {
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    frictionConeResidual(mu_, s.f[i], data.residual.template segment<5>(5*i));
    data.slack.template segment<5>(5*i)
        = - data.residual.template segment<5>(5*i);
  }
  setSlackAndDualPositive(data);
}


void LinearizedFrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dtau * Jac_.transpose() * data.dual.template segment<5>(5*i);
      dimf_stack += 3;
    }
  }
}


void LinearizedFrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx = 5*i;
      data.r[i].array() 
          = (data.dual.template segment<5>(idx).array()
              *data.residual.template segment<5>(idx).array()
              -data.duality.template segment<5>(idx).array())
              / data.slack.template segment<5>(idx).array();
      kkt_residual.lf().template segment<3>(dimf_stack).noalias()
          += dtau * Jac_.transpose() * data.r[i];
      data.r[i].array() = data.dual.template segment<5>(idx).array() 
                          / data.slack.template segment<5>(idx).array();
      data.J[i].noalias() = data.r[i].asDiagonal() * Jac_;
      kkt_matrix.Qff().template block<3, 3>(dimf_stack, dimf_stack).noalias()
          += dtau * Jac_.transpose() * data.J[i];
      dimf_stack += 3;
    }
  }
}


void LinearizedFrictionCone::computeSlackAndDualDirection(
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
      const int idx = 5*i;
      data.dslack.template segment<5>(idx).noalias()
          = - Jac_ * d.df().segment<3>(dimf_stack) 
              - data.residual.template segment<5>(idx);
      for (int j=0; j<5; ++j) {
        data.ddual.coeffRef(idx+j) = computeDualDirection(data.slack.coeff(idx+j), 
                                                          data.dual.coeff(idx+j), 
                                                          data.dslack.coeff(idx+j), 
                                                          data.duality.coeff(idx+j));
      }
      dimf_stack += 3;
    }
  }
}


void LinearizedFrictionCone::computePrimalAndDualResidual(
    Robot& robot, ConstraintComponentData& data, const SplitSolution& s) const {
  data.residual.setZero();
  data.duality.setZero();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      const int idx = 5*i;
      frictionConeResidual(mu_, s.f[i], data.residual.template segment<5>(idx));
      data.residual.template segment<5>(idx).noalias()
          += data.slack.template segment<5>(idx);
      for (int j=0; j<5; ++j) {
        data.duality.coeffRef(idx+j) = computeDuality(data.slack.coeff(idx+j), 
                                                      data.dual.coeff(idx+j));
      }
    }
  }
}


int LinearizedFrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp