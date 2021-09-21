#ifndef IDOCP_SWITCHING_TIME_COST_FUNCTION_HXX_
#define IDOCP_SWITCHING_TIME_COST_FUNCTION_HXX_

#include "idocp/hybrid/switching_time_cost_function.hpp"

namespace idocp {

inline SwitchingTimeCostFunction::SwitchingTimeCostFunction()
  : costs_(),
    ts_(),
    hts_(),
    Qts_() {
}


inline SwitchingTimeCostFunction::~SwitchingTimeCostFunction() {
}


inline void SwitchingTimeCostFunction::push_back(
    const SwitchingTimeCostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


inline void SwitchingTimeCostFunction::clear() {
  costs_.clear();
}


inline double SwitchingTimeCostFunction::computeCost(
    const HybridTimeDiscretization& discretization) {
  if (!costs_.empty()) {
    const int num_events = discretization.N_impulse() + discretization.N_lift();
    setNumSwitches(num_events);
    setSwitchingTimes(discretization);
    const double t0 = discretization.t(0);
    const double tf = discretization.t(discretization.N());
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->computeCost(t0, tf, ts_);
    }
    return cost;
  }
  else {
    return 0;
  }
}


inline double SwitchingTimeCostFunction::linearizeCost(
    const HybridTimeDiscretization& discretization, KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = discretization.N_impulse() + discretization.N_lift();
    setNumSwitches(num_events);
    setSwitchingTimes(discretization);
    const double t0 = discretization.t(0);
    const double tf = discretization.t(discretization.N());
    double cost = 0;
    hts_.setZero();
    for (auto& e : costs_) {
      cost += e->computeCost(t0, tf, ts_);
      e->computeCostDerivatives(t0, tf, ts_, hts_);
    }
    return cost;
  }
  else {
    return 0;
  }
}


inline double SwitchingTimeCostFunction::quadratizeStageCost(
    const HybridTimeDiscretization& discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = discretization.N_impulse() + discretization.N_lift();
    setNumSwitches(num_events);
    setSwitchingTimes(discretization);
    const double t0 = discretization.t(0);
    const double tf = discretization.t(discretization.N());
    double cost = 0;
    hts_.setZero();
    Qts_.setZero();
    for (auto& e : costs_) {
      cost += e->computeCost(t0, tf, ts_);
      e->computeCostDerivatives(t0, tf, ts_, hts_);
      e->computeCostHessian(t0, tf, ts_, Qts_);
    }
    return cost;
  }
  else {
    return 0;
  }
}


inline void SwitchingTimeCostFunction::setNumSwitches(const int num_switches) {
  if (ts_.size() != num_switches) {
    ts_.resize(num_switches);
  }
  if (hts_.size() != num_switches) {
    hts_.resize(num_switches);
  }
  if (Qts_.cols() != num_switches) {
    Qts_.resize(num_switches, num_switches);
  }
}


inline void SwitchingTimeCostFunction::setSwitchingTimes(
    const HybridTimeDiscretization& discretization) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    ts_.coeffRef(0) = discretization.t_impulse(0);
    ++impulse_index;
  }
  else {
    ts_.coeffRef(0) = discretization.t_lift(0);
    ++lift_index;
  }
  for (int event_index=1; event_index<num_events; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      ts_.coeffRef(event_index) = discretization.t_impulse(impulse_index);
      ++impulse_index;
    }
    else {
      ts_.coeffRef(event_index) = discretization.t_lift(lift_index);
      ++lift_index;
    }
  }
}


inline void SwitchingTimeCostFunction::setKKT(
    const HybridTimeDiscretization& discretization, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    kkt_residual.aux[0].h -= hts.coeff(0);
    kkt_matrix.aux[0].Qtt += Qts.coeff(0, 0);
    ++impulse_index;
  }
  else {
    kkt_residual.lift[0].h -= hts.coeff(0);
    kkt_matrix.lift[0].Qtt += Qts.coeff(0, 0);
    ++lift_index;
  }
  for (int event_index=1; event_index<num_events; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      kkt_residual.aux[impulse_index].h -= hts.coeff(event_index);
      kkt_matrix.aux[impulse_index].Qtt += Qts.coeff(event_index, event_index);
      ++impulse_index;
    }
    else {
      kkt_residual.lift[lift_index].h -= hts.coeff(event_index);
      kkt_matrix.lift[lift_index].Qtt += Qts.coeff(event_index, event_index);
      ++lift_index;
    }
  }
}


inline void SwitchingTimeCostFunction::setKKT(
    const HybridTimeDiscretization& discretization, KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  int impulse_index = 0;
  int lift_index = 0;
  if (discretization.eventType(0) == DiscreteEventType::Impulse) {
    kkt_residual.aux[0].h -= hts.coeff(0);
    ++impulse_index;
  }
  else {
    kkt_residual.lift[0].h -= hts.coeff(0);
    ++lift_index;
  }
  for (int event_index=1; event_index<num_events; ++event_index) {
    assert(event_index == impulse_index+lift_index);
    const auto next_event_type = discretization.eventType(event_index+1);
    if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
      kkt_residual.aux[impulse_index].h -= hts.coeff(event_index);
      ++impulse_index;
    }
    else {
      kkt_residual.lift[lift_index].h -= hts.coeff(event_index);
      ++lift_index;
    }
  }
}

} // namespace idocp

#endif // IDOCP_SWITCHING_TIME_COST_FUNCTION_HXX_ 