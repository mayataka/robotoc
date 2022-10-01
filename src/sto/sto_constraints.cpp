#include "robotoc/sto/sto_constraints.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

STOConstraints::STOConstraints(const int reserved_num_switches, 
                               const double min_dt, 
                               const double barrier_param, 
                               const double fraction_to_boundary_rule) 
  : STOConstraints(std::vector<double>(reserved_num_switches+1, min_dt), 
                   barrier_param, fraction_to_boundary_rule) {
}


STOConstraints::STOConstraints(const std::vector<double>& min_dt, 
                               const double barrier_param, 
                               const double fraction_to_boundary_rule) 
  : dtlb_(min_dt.size(), DwellTimeLowerBound(barrier_param, 
                                             fraction_to_boundary_rule)),
    min_dt_(min_dt), 
    eps_(std::sqrt(std::numeric_limits<double>::epsilon())),
    barrier_(barrier_param), 
    fraction_to_boundary_rule_(fraction_to_boundary_rule),
    reserved_num_switches_(min_dt.size()-1),
    num_switches_(0),
    primal_step_size_(Eigen::VectorXd::Zero(min_dt.size())), 
    dual_step_size_(Eigen::VectorXd::Zero(min_dt.size())) {
  for (const auto e : min_dt) {
    if (e < 0.) {
      throw std::out_of_range(
          "[STOConstraints] invalid argment: 'min_dt' must be non-negative!");
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
  minimum_dwell_times_.resize(min_dt.size());
  for (int i=0; i<min_dt_.size(); ++i) {
    minimum_dwell_times_[i] = min_dt[i];
  }
}


STOConstraints::STOConstraints()
  : dtlb_(),
    min_dt_(0), 
    num_switches_(0),
    primal_step_size_(), 
    dual_step_size_() {
}


STOConstraints::~STOConstraints() {
}


void STOConstraints::setSlack(const TimeDiscretization& time_discretization) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  reserve(num_events);
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
    dtlb_[0].setSlack(min_dt_[0], time_discretization.t0(), time_discretization.impulseTime(0));
  }
  else {
    dtlb_[0].setSlack(min_dt_[0], time_discretization.t0(), time_discretization.liftTime(0));
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = time_discretization.eventType(event_index+1);
    if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      time_discretization.impulseTime(impulse_index),  
                                      time_discretization.impulseTime(impulse_index+1));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1], 
                                      time_discretization.impulseTime(impulse_index),  
                                      time_discretization.liftTime(lift_index));
      }
      ++impulse_index;
    }
    else {
      if (next_event_type == DiscreteEventType::Impulse) {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      time_discretization.liftTime(lift_index),  
                                      time_discretization.impulseTime(impulse_index));
      }
      else {
        dtlb_[event_index+1].setSlack(min_dt_[event_index+1],
                                      time_discretization.liftTime(lift_index),  
                                      time_discretization.liftTime(lift_index+1));
      }
      ++lift_index;
    }
  }
  if (time_discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    dtlb_[num_events].setSlack(min_dt_[num_events], 
                               time_discretization.impulseTime(impulse_index), 
                               time_discretization.tf());
  }
  else {
    dtlb_[num_events].setSlack(min_dt_[num_events], 
                               time_discretization.liftTime(lift_index), 
                               time_discretization.tf());
  }
}


void STOConstraints::evalConstraint(const TimeDiscretization& time_discretization) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  reserve(num_events);
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], time_discretization.t0(), 
                              time_discretization.impulseTime(0));
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], time_discretization.t0(), 
                              time_discretization.liftTime(0));
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = time_discretization.eventType(event_index+1);
    if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.impulseTime(impulse_index), 
              time_discretization.impulseTime(impulse_index+1));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.impulseTime(impulse_index), 
              time_discretization.liftTime(lift_index));
        }
      }
      ++impulse_index;
    }
    else {
      if (time_discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.liftTime(lift_index), 
              time_discretization.impulseTime(impulse_index));
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.liftTime(lift_index), 
              time_discretization.liftTime(lift_index+1));
        }
      }
      ++lift_index;
    }
  }
  if (time_discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       time_discretization.impulseTime(impulse_index), 
                                       time_discretization.tf());
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       time_discretization.liftTime(lift_index), 
                                       time_discretization.tf());
    }
  }
}


void STOConstraints::linearizeConstraints(
    const TimeDiscretization& time_discretization, KKTResidual& kkt_residual) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  reserve(num_events);
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], time_discretization.t0(), 
                              time_discretization.impulseTime(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.aux[0]);
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(0)) {
      dtlb_[0].evalConstraint(min_dt_[0], time_discretization.t0(), 
                              time_discretization.liftTime(0));
      dtlb_[0].evalDerivatives_lb(kkt_residual.lift[0]);
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = time_discretization.eventType(event_index+1);
    if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.impulseTime(impulse_index), 
              time_discretization.impulseTime(impulse_index+1));
          dtlb_[event_index+1].evalDerivatives_lub( 
              kkt_residual.aux[impulse_index], 
              kkt_residual.aux[impulse_index+1]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.impulseTime(impulse_index), 
              time_discretization.liftTime(lift_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.aux[impulse_index], kkt_residual.lift[lift_index]);
        }
      }
      ++impulse_index;
    }
    else {
      if (time_discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].evalConstraint(
               min_dt_[event_index+1], time_discretization.liftTime(lift_index), 
              time_discretization.impulseTime(impulse_index));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.aux[impulse_index]);
        }
        else {
          dtlb_[event_index+1].evalConstraint(
              min_dt_[event_index+1], time_discretization.liftTime(lift_index), 
              time_discretization.liftTime(lift_index+1));
          dtlb_[event_index+1].evalDerivatives_lub(
              kkt_residual.lift[lift_index], kkt_residual.lift[lift_index+1]);
        }
      }
      ++lift_index;
    }
  }
  if (time_discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       time_discretization.impulseTime(impulse_index), 
                                       time_discretization.tf());
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].evalConstraint(min_dt_[num_events], 
                                       time_discretization.liftTime(lift_index), 
                                       time_discretization.tf());
      dtlb_[num_events].evalDerivatives_ub(kkt_residual.lift[lift_index]);
    }
  }
}


void STOConstraints::condenseSlackAndDual(
    const TimeDiscretization& time_discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  reserve(num_events);
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].condenseSlackAndDual_lb(kkt_matrix.aux[0], kkt_residual.aux[0]);
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(0)) {
      dtlb_[0].condenseSlackAndDual_lb(kkt_matrix.lift[0], kkt_residual.lift[0]);
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = time_discretization.eventType(event_index+1);
    if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
              kkt_matrix.aux[impulse_index+1], kkt_residual.aux[impulse_index+1]);
        }
        else {
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
        }
      }
      ++impulse_index;
    }
    else {
      if (time_discretization.isSTOEnabledLift(lift_index)) {
        if (next_event_type == DiscreteEventType::Impulse) {
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
              kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
        }
        else {
          dtlb_[event_index+1].condenseSlackAndDual_lub(
              kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index],
              kkt_matrix.lift[lift_index+1], kkt_residual.lift[lift_index+1]);
        }
      }
      ++lift_index;
    }
  }
  if (time_discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index]);
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].condenseSlackAndDual_ub(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index]);
    }
  }
}


void STOConstraints::expandSlackAndDual(
    const TimeDiscretization& time_discretization, const Direction& d) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  reserve(num_events);
  num_switches_ = num_events;
  assert(num_events+1 <= dtlb_.size());
  if (num_events <= 0) {
    return;
  } 
  primal_step_size_.fill(1.0);
  dual_step_size_.fill(1.0);
  if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.aux[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
  }
  else {
    if (time_discretization.isSTOEnabledLift(0)) {
      dtlb_[0].expandSlackAndDual_lb(d.lift[0].dts);
      primal_step_size_.coeffRef(0) = dtlb_[0].maxSlackStepSize();
      dual_step_size_.coeffRef(0) = dtlb_[0].maxDualStepSize();
    }
  }
  int impulse_index = 0;
  int lift_index = 0;
  for (int event_index=0; event_index<num_events-1; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = time_discretization.eventType(event_index+1);
    if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
        dtlb_[event_index+1].expandSlackAndDual_lub(d.aux[impulse_index].dts,
                                                    d.aux[impulse_index].dts_next);
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      ++impulse_index;
    }
    else {
      if (time_discretization.isSTOEnabledLift(lift_index)) {
        dtlb_[event_index+1].expandSlackAndDual_lub(d.lift[lift_index].dts,
                                                    d.lift[lift_index].dts_next);
        primal_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxSlackStepSize();
        dual_step_size_.coeffRef(event_index+1) 
            = dtlb_[event_index+1].maxDualStepSize();
      }
      ++lift_index;
    }
  }
  if (time_discretization.eventType(num_events-1) == DiscreteEventType::Impulse) {
    if (time_discretization.isSTOEnabledImpulse(impulse_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(d.aux[impulse_index].dts);
      primal_step_size_.coeffRef(num_events) = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) = dtlb_[num_events].maxDualStepSize();
    } 
  }
  else {
    if (time_discretization.isSTOEnabledLift(lift_index)) {
      dtlb_[num_events].expandSlackAndDual_ub(d.lift[lift_index].dts);
      primal_step_size_.coeffRef(num_events) = dtlb_[num_events].maxSlackStepSize();
      dual_step_size_.coeffRef(num_events) = dtlb_[num_events].maxDualStepSize();
    }
  }
}


double STOConstraints::maxPrimalStepSize() const {
  if (num_switches_ > 0) {
    return primal_step_size_.head(num_switches_+1).minCoeff();
  }
  else {
    return 1.0;
  }
}


double STOConstraints::maxDualStepSize() const {
  if (num_switches_ > 0) {
    return dual_step_size_.head(num_switches_+1).minCoeff();
  }
  else {
    return 1.0;
  }
}


void STOConstraints::updateSlack(const double step_size) {
  for (int i=0; i<num_switches_+1; ++i) {
    dtlb_[i].updateSlack(step_size);
  }
}


void STOConstraints::updateDual(const double step_size) {
  for (int i=0; i<num_switches_+1; ++i) {
    dtlb_[i].updateDual(step_size);
  }
}


double STOConstraints::KKTError() const {
  double err = 0;
  for (int i=0; i<num_switches_+1; ++i) {
    err += dtlb_[i].KKTError();
  }
  return err;
}


void STOConstraints::setMinimumDwellTimes(const double min_dt) {
  if (min_dt < 0) {
    throw std::out_of_range(
        "[STOConstraints] invalid argment: 'min_dt' must be non-negative!");
  }
  for (auto& e : min_dt_) {
    e = min_dt;
  }
}


void STOConstraints::setMinimumDwellTimes(
    const std::vector<double>& min_dt) {
  for (const auto e : min_dt) {
    if (e < 0) {
      throw std::out_of_range(
          "[STOConstraints] invalid argment: 'min_dt' must be non-negative!");
    }
  }
  min_dt_ = min_dt;
  while (min_dt_.size() < (reserved_num_switches_+1)) {
    min_dt_.push_back(eps_);
  }
  reserve(min_dt_.size()-1);
}


const std::vector<double>& STOConstraints::getMinimumDwellTimes() const {
  return min_dt_;
}


void STOConstraints::setBarrierParam(const double barrier_param) {
  assert(barrier_param > 0.0);
  for (auto& e : dtlb_) {
    e.setBarrierParam(barrier_param);
  }
  barrier_ = barrier_param;
}


void STOConstraints::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  assert(fraction_to_boundary_rule > 0.0);
  assert(fraction_to_boundary_rule < 1.0);
  for (auto& e : dtlb_) {
    e.setFractionToBoundaryRule(fraction_to_boundary_rule);
  }
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


double STOConstraints::getBarrierParam() const {
  return barrier_;
}


double STOConstraints::getFractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


void STOConstraints::reserve(const int reserved_num_switches) { 
  if (reserved_num_switches_ < reserved_num_switches) {
    while (dtlb_.size() < (reserved_num_switches+1)) {
      dtlb_.emplace_back(barrier_, fraction_to_boundary_rule_);
    }
    while (min_dt_.size() < (reserved_num_switches+1)) {
      min_dt_.push_back(eps_);
    }
    if (primal_step_size_.size() < reserved_num_switches+1) {
      primal_step_size_.conservativeResize(reserved_num_switches+1);
    }
    if (dual_step_size_.size() < reserved_num_switches+1) {
      dual_step_size_.conservativeResize(reserved_num_switches+1);
    }
    reserved_num_switches_ = reserved_num_switches;
  }
}


int STOConstraints::reservedNumSwitches() const {
  return reserved_num_switches_;
}



ConstraintComponentData STOConstraints::createConstraintsData(
    const TimeDiscretization& time_discretization) const {
  const int N = time_discretization.N_grids();
  const int num_discrete_events = time_discretization.grid(N).contact_phase 
                                    - time_discretization.grid(0).contact_phase;
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
  data.r[1].resize(num_phases);
  for (int i=0; i<num_phases; ++i) {
    data.r[1].coeffRef(i) = min_dt_[i];
  }
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
  // pdipm::setSlackAndDualPositive(barrier_, data);
  for (int i=0; i<data.slack.size(); ++i) {
    if (data.slack.coeff(i) < 0.0) {
      data.slack.coeffRef(i) = std::sqrt(barrier_);
    }
  }
  data.dual.array() = barrier_ / data.slack.array();
}


void STOConstraints::evalConstraint(
    const TimeDiscretization& time_discretization, 
    ConstraintComponentData& data) const {
  computeDwellTimes(time_discretization, data.r[0]);
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
  return pdipm::fractionToBoundaryDual(fraction_to_boundary_rule_, data);
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
  const int N = time_discretization.N_grids();
  const int num_discrete_events = time_discretization.grid(N).contact_phase 
                                    - time_discretization.grid(0).contact_phase;
  dwell_times.resize(num_discrete_events+1);
  dwell_times.setZero();
  int event_index = 0;
  double prev_event_time = time_discretization.grid(0).t;
  for (int i=0; i<N; ++i) {
    if (time_discretization.grid(i).type == GridType::Impulse
        || time_discretization.grid(i).type == GridType::Lift) {
      dwell_times.coeffRef(event_index) = time_discretization.grid(i).t - prev_event_time;
      prev_event_time = time_discretization.grid(i).t;
      ++event_index;
    }
  }
  dwell_times.coeffRef(num_discrete_events) 
      = time_discretization.grid(N).t - prev_event_time;
}

} // namespace robotoc
