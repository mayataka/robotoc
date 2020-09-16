#include "idocp/constraints/friction_cone.hpp"

#include <exception>
#include <iostream>
#include <assert.h>


namespace idocp {

FrictionCone::FrictionCone(
    const Robot& robot, const int contact_index, const double mu,
    const double barrier, const double fraction_to_boundary_rate)
  : ConstraintComponentBase(barrier, fraction_to_boundary_rate),
    contact_index_(contact_index),
    dimc_(robot.jointVelocityLimit().size()),
    mu_(mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range("invalid argment: mu must be positive");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


FrictionCone::FrictionCone()
  : ConstraintComponentBase(),
    contact_index_(0),
    dimc_(0),
    mu_(0) {
}


FrictionCone::~FrictionCone() {
}


bool FrictionCone::useKinematics() const {
  return false;
}


bool FrictionCone::isFeasible(Robot& robot, ConstraintComponentData& data, 
                              const SplitSolution& s) const {
  s.f.
  if (s.v.tail(dimc_).coeff(i) > vmax_.coeff(i)) {
    return false;
  }
  return true;
}


void FrictionCone::setSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau > 0);
  data.slack = dtau * (vmax_-s.v.tail(dimc_));
  setSlackAndDualPositive(data.slack, data.dual);
}


void FrictionCone::augmentDualResidual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    KKTResidual& kkt_residual) const {
  kkt_residual.lv().tail(dimc_).noalias() += dtau * data.dual;
}


void FrictionCone::condenseSlackAndDual(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) const {
  kkt_matrix.Qvv().diagonal().tail(dimc_).array()
      += dtau * dtau * data.dual.array() / data.slack.array();
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  kkt_residual.lv().tail(dimc_).array() 
      += dtau * (data.dual.array()*data.residual.array()-data.duality.array()) 
              / data.slack.array();
}


void FrictionCone::computeSlackAndDualDirection(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitDirection& d) const {
  data.dslack = - dtau * d.dv().tail(dimc_) - data.residual;
  computeDualDirection(data.slack, data.dual, data.dslack, data.duality, 
                       data.ddual);
}


double FrictionCone::residualL1Nrom(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  return data.residual.lpNorm<1>();
}


double FrictionCone::squaredKKTErrorNorm(
    Robot& robot, ConstraintComponentData& data, const double dtau, 
    const SplitSolution& s) const {
  data.residual = dtau * (s.v.tail(dimc_)-vmax_) + data.slack;
  computeDuality(data.slack, data.dual, data.duality);
  double error = 0;
  error += data.residual.squaredNorm();
  error += data.duality.squaredNorm();
  return error;
}


int FrictionCone::dimc() const {
  return dimc_;
}

} // namespace idocp