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

} // namespace robotoc
