#include "robotoc/sto/sto_constraints.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

STOConstraints::STOConstraints(const int num_switches, 
                               const double minimum_dwell_time, 
                               const double barrier_param, 
                               const double fraction_to_boundary_rule) 
  : STOConstraints(std::vector<double>(num_switches+1, minimum_dwell_time), 
                   barrier_param, fraction_to_boundary_rule) {
}


STOConstraints::STOConstraints(const std::vector<double>& minimum_dwell_times, 
                               const double barrier_param, 
                               const double fraction_to_boundary_rule) 
  : STOConstraints(Eigen::Map<const Eigen::VectorXd>(minimum_dwell_times.data(), minimum_dwell_times.size()), 
                   barrier_param, fraction_to_boundary_rule) {
}


STOConstraints::STOConstraints(const Eigen::VectorXd& minimum_dwell_times, 
                               const double barrier_param, 
                               const double fraction_to_boundary_rule) 
  : barrier_(barrier_param), 
    fraction_to_boundary_rule_(fraction_to_boundary_rule),
    num_switches_(minimum_dwell_times.size()-1),
    primal_step_size_(Eigen::VectorXd::Zero(minimum_dwell_times.size())), 
    dual_step_size_(Eigen::VectorXd::Zero(minimum_dwell_times.size())),
    minimum_dwell_times_(minimum_dwell_times) {
  for (int i=0; i<minimum_dwell_times.size(); ++i) {
    if (minimum_dwell_times.coeff(i) < 0.) {
      throw std::out_of_range(
          "[STOConstraints] invalid argment: 'minimum_dwell_times' must be non-negative!");
    }
  }
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[STOConstraints] invalid argment: 'barrier_param' must be positive!");
  }
  if (fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[STOConstraints] invalid argment: 'fraction_to_boundary_rule' must be positive!");
  }
  if (fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[STOConstraints] invalid argment: 'fraction_to_boundary_rule' must be less than 1!");
  }
}


STOConstraints::STOConstraints()
  : barrier_(1.0e-03), 
    fraction_to_boundary_rule_(0.995),
    num_switches_(0),
    primal_step_size_(), 
    dual_step_size_(),
    minimum_dwell_times_() {
}


void STOConstraints::setMinimumDwellTimes(const double minimum_dwell_time) {
  setMinimumDwellTimes(std::vector<double>(num_switches_+1, minimum_dwell_time));
}


void STOConstraints::setMinimumDwellTimes(
    const std::vector<double>& minimum_dwell_times) {
  setMinimumDwellTimes(Eigen::Map<const Eigen::VectorXd>(minimum_dwell_times.data(), 
                                                         minimum_dwell_times.size()));
}


void STOConstraints::setMinimumDwellTimes(
    const Eigen::VectorXd& minimum_dwell_times) {
  for (int i=0; i<minimum_dwell_times.size(); ++i) {
    if (minimum_dwell_times.coeff(i) < 0.) {
      throw std::out_of_range(
          "[STOConstraints] invalid argment: 'minimum_dwell_times' must be non-negative!");
    }
  }
  minimum_dwell_times_ = minimum_dwell_times;
  num_switches_ = minimum_dwell_times.size() - 1;
}


const Eigen::VectorXd& STOConstraints::getMinimumDwellTimes() const {
  return minimum_dwell_times_;
}


void STOConstraints::setBarrierParam(const double barrier_param) {
  assert(barrier_param > 0.0);
  barrier_ = barrier_param;
}


void STOConstraints::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  assert(fraction_to_boundary_rule > 0.0);
  assert(fraction_to_boundary_rule < 1.0);
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


double STOConstraints::getBarrierParam() const {
  return barrier_;
}


double STOConstraints::getFractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


ConstraintComponentData STOConstraints::createConstraintsData(
    const TimeDiscretization& time_discretization) const {
  const int N = time_discretization.size() - 1;
  const int num_discrete_events = time_discretization.grid(N).phase 
                                    - time_discretization.grid(0).phase;
  auto data = ConstraintComponentData(num_discrete_events+1, barrier_);
  data.r.push_back(Eigen::VectorXd::Zero(num_discrete_events+1));
  data.r.push_back(Eigen::VectorXd::Zero(num_discrete_events+1));
  data.J.push_back(Eigen::MatrixXd::Zero(num_discrete_events+1, num_discrete_events));
  data.J.push_back(Eigen::MatrixXd::Zero(num_discrete_events+1, num_discrete_events));
  return data;
}


bool STOConstraints::isFeasible(const TimeDiscretization& time_discretization,
                                ConstraintComponentData& data) const {
  computeDwellTimes(time_discretization, data.r[0]);
  if (data.r[0].size() != minimum_dwell_times_.size()) {
    throw std::runtime_error(
        "[STOConstraints] : invalid size of minimum_dwell_times_ is detected! It should be " + std::to_string(data.r[0].size()));
  }
  data.residual = minimum_dwell_times_ - data.r[0];
  return (data.residual.minCoeff() > 0.0);
}


void STOConstraints::setSlackAndDual(
    const TimeDiscretization& time_discretization, 
    ConstraintComponentData& data) const {
  const int num_phases = minimum_dwell_times_.size();
  if (num_phases <= 1) return;

  data.r[0].resize(num_phases);
  computeDwellTimes(time_discretization, data.r[0]);
  if (data.r[0].size() != minimum_dwell_times_.size()) {
    throw std::runtime_error(
        "[STOConstraints] : invalid size of minimum_dwell_times_ is detected! It should be " + std::to_string(data.r[0].size()));
  }
  data.r[1].resize(num_phases);
  data.r[1] = minimum_dwell_times_;
  data.J[0].resize(num_phases, num_phases-1);
  data.J[1].resize(num_phases, num_phases-1);
  auto& J = data.J[0];
  J.setZero();
  J.coeffRef(0, 0) = -1.0;
  for (int i=0; i<num_phases-2; ++i) {
    J.coeffRef(i+1, i)   =  1.0;
    J.coeffRef(i+1, i+1) = -1.0;
  }
  J.coeffRef(num_phases-1, num_phases-2) = 1.0;

  data.resize(num_phases);
  data.slack = - (minimum_dwell_times_ - data.r[0]);
  pdipm::setSlackAndDualPositive(barrier_, data);
}


void STOConstraints::evalConstraint(
    const TimeDiscretization& time_discretization, 
    ConstraintComponentData& data) const {
  computeDwellTimes(time_discretization, data.r[0]);
  if (data.r[0].size() != minimum_dwell_times_.size()) {
    throw std::runtime_error(
        "[STOConstraints] : invalid size of minimum_dwell_times_ is detected! It should be " + std::to_string(data.r[0].size()));
  }
  data.residual = minimum_dwell_times_ - data.r[0] + data.slack;
  pdipm::computeComplementarySlackness(barrier_, data);
  data.log_barrier = pdipm::logBarrier(barrier_, data.slack);
}


void STOConstraints::linearizeConstraints(
    const TimeDiscretization& time_discretization,  
    ConstraintComponentData& data, Eigen::VectorXd& lt) const {
  evalConstraint(time_discretization, data);
  lt.noalias() += data.J[0].transpose() * data.dual;
}


void STOConstraints::condenseSlackAndDual(ConstraintComponentData& data, 
                                          Eigen::VectorXd& lt,
                                          Eigen::MatrixXd& Qtt) const {
  pdipm::computeCondensingCoeffcient(data);
  lt.noalias() += data.J[0].transpose() * data.cond;
  data.cond.array() = data.dual.array() / data.slack.array();
  data.J[1].noalias() = data.cond.asDiagonal() * data.J[0];
  Qtt.noalias() += data.J[0].transpose() * data.J[1];
}


void STOConstraints::expandSlackAndDual(ConstraintComponentData& data, 
                                        Eigen::VectorXd& dts) const {
  data.dslack.noalias() = - data.J[0] * dts - data.residual;
  pdipm::computeDualDirection(data);
}


double STOConstraints::maxSlackStepSize(const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


double STOConstraints::maxDualStepSize(const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundaryDual(fraction_to_boundary_rule_, data);
}


void STOConstraints::updateSlack(ConstraintComponentData& data, 
                                 const double step_size) const {
  assert(step_size > 0);
  assert(step_size <= 1.0);
  data.slack.noalias() += step_size * data.dslack;
}


void STOConstraints::updateDual(ConstraintComponentData& data, 
                                const double step_size) const {
  assert(step_size > 0);
  assert(step_size <= 1.0);
  data.dual.noalias() += step_size * data.ddual;
}


void STOConstraints::computeDwellTimes(const TimeDiscretization& time_discretization,
                                       Eigen::VectorXd& dwell_times) {
  const int N = time_discretization.size() - 1;
  const int num_discrete_events = time_discretization.grid(N).phase 
                                    - time_discretization.grid(0).phase;
  dwell_times.resize(num_discrete_events+1);
  dwell_times.setZero();
  int event_index = 0;
  double prev_event_time = time_discretization.grid(0).t;
  for (int i=0; i<N; ++i) {
    if (time_discretization[i].type == GridType::Impact
        || time_discretization[i].type == GridType::Lift) {
      dwell_times.coeffRef(event_index) = time_discretization[i].t - prev_event_time;
      prev_event_time = time_discretization[i].t;
      ++event_index;
    }
  }
  dwell_times.coeffRef(num_discrete_events) 
      = time_discretization.grid(N).t - prev_event_time;
}

} // namespace robotoc
