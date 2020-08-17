#include "idocp/constraints/constraint_component_base.hpp"


namespace idocp {

ConstraintComponentBase::ConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rate) 
  : barrier_(barrier),
    fraction_to_boundary_rate_(fraction_to_boundary_rate) {
}


ConstraintComponentBase::ConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rate_(0) {
}


void ConstraintComponentBase::setBarrier(const double barrier) {
  barrier_ = barrier;
}


void ConstraintComponentBase::setFractionToBoundaryRate(
    const double fraction_to_boundary_rate) {
  fraction_to_boundary_rate_ = fraction_to_boundary_rate;
}


double ConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return fractionToBoundary(data.slack, data.dslack);
}


double ConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return fractionToBoundary(data.dual, data.ddual);
}


void ConstraintComponentBase::updateSlack(ConstraintComponentData& data, 
                                          const double step_size) const {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


void ConstraintComponentBase::updateDual(ConstraintComponentData& data, 
                                         const double step_size) const {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data) const {
  const double cost = - barrier_ * data.slack.array().log().sum();
  return cost;
}


double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  const double cost 
      = - barrier_ * (data.slack+step_size*data.dslack).array().log().sum();
  return cost;
}


void ConstraintComponentBase::setSlackAndDualPositive(
    Eigen::VectorXd& slack, Eigen::VectorXd& dual) const {
  for (int i=0; i<dimc(); ++i) {
    while (slack.coeff(i) < barrier_) {
      slack.coeffRef(i) += barrier_;
    }
  }
  dual.array() = barrier_ / slack.array();
}


void ConstraintComponentBase::computeDualityResidual(
    const Eigen::VectorXd& slack, const Eigen::VectorXd& dual, 
    Eigen::VectorXd& duality) const {
  duality.array() = slack.array() * dual.array() - barrier_;
}


void ConstraintComponentBase::computeDualDirection(
    const Eigen::VectorXd& slack, const Eigen::VectorXd& dslack, 
    const Eigen::VectorXd& dual, const Eigen::VectorXd& duality, 
    Eigen::VectorXd& ddual) const {
  ddual.array() = - (dual.array()*dslack.array()+duality.array())/slack.array();
}


double ConstraintComponentBase::fractionToBoundary(
    const Eigen::VectorXd& vec, const Eigen::VectorXd& dvec) const {
  double min_fraction_to_boundary = 1;
  for (int i=0; i<dimc(); ++i) {
    const double fraction_to_boundary 
        = - fraction_to_boundary_rate_ * (vec.coeff(i)/dvec.coeff(i));
    if (fraction_to_boundary > 0 && fraction_to_boundary < 1) {
      if (fraction_to_boundary < min_fraction_to_boundary) {
        min_fraction_to_boundary = fraction_to_boundary;
      }
    }
  }
  return min_fraction_to_boundary;
}

} // namespace idocp
