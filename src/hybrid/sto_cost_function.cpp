#include "robotoc/hybrid/sto_cost_function.hpp"


namespace robotoc {

STOCostFunction::STOCostFunction()
  : costs_(),
    lts_(),
    Qts_() {
}


STOCostFunction::~STOCostFunction() {
}


void STOCostFunction::push_back(const STOCostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


void STOCostFunction::clear() {
  costs_.clear();
}


double STOCostFunction::evalCost(const TimeDiscretization& time_discretization) {
  if (!costs_.empty()) {
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->evalCost(time_discretization);
    }
    return cost;
  }
  else {
    return 0.0;
  }
}


double STOCostFunction::linearizeCost(const TimeDiscretization& time_discretization, 
                                      KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
    lts_.resize(num_events);
    lts_.setZero();
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->evalCost(time_discretization);
      e->evalCostDerivatives(time_discretization, lts_);
    }
    setToKKT(time_discretization, kkt_residual);
    return cost;
  }
  else {
    return 0.0;
  }
}


double STOCostFunction::quadratizeCost(const TimeDiscretization& time_discretization, 
                                       KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
    if (lts_.size() != num_events) {
      lts_.resize(num_events);
    }
    if ((Qts_.cols() != num_events) || (Qts_.rows() != num_events)) {
      Qts_.resize(num_events, num_events);
    }
    lts_.setZero();
    Qts_.setZero();
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->evalCost(time_discretization);
      e->evalCostDerivatives(time_discretization, lts_);
      e->evalCostHessian(time_discretization, Qts_);
    }
    setToKKT(time_discretization, kkt_matrix, kkt_residual);
    return cost;
  }
  else {
    return 0.0;
  }
}


void STOCostFunction::setToKKT(const TimeDiscretization& time_discretization, 
                               KKTResidual& kkt_residual) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
      kkt_residual.aux[0].h -= lts_.coeff(0);
      ++impulse_index;
    }
    else {
      kkt_residual.lift[0].h -= lts_.coeff(0);
      ++lift_index;
    }
    for (int event_index=1; event_index<num_events; ++event_index) {
      assert(event_index == impulse_index+lift_index);
      if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
        kkt_residual.aux[impulse_index].h -= lts_.coeff(event_index);
        ++impulse_index;
      }
      else {
        kkt_residual.lift[lift_index].h -= lts_.coeff(event_index);
        ++lift_index;
      }
    }
  }
}


void STOCostFunction::setToKKT(const TimeDiscretization& time_discretization, 
                               KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int num_events = time_discretization.N_impulse() + time_discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (time_discretization.eventType(0) == DiscreteEventType::Impulse) {
      kkt_residual.aux[0].h -= lts_.coeff(0);
      kkt_matrix.aux[0].Qtt += Qts_.coeff(0, 0);
      kkt_matrix.aux[0].Qtt_prev += Qts_.coeff(0, 1);
      ++impulse_index;
    }
    else {
      kkt_residual.lift[0].h -= lts_.coeff(0);
      kkt_matrix.lift[0].Qtt += Qts_.coeff(0, 0);
      kkt_matrix.lift[0].Qtt_prev += Qts_.coeff(0, 1);
      ++lift_index;
    }
    for (int event_index=1; event_index<num_events; ++event_index) {
      assert(event_index == impulse_index+lift_index);
      if (time_discretization.eventType(event_index) == DiscreteEventType::Impulse) {
        kkt_residual.aux[impulse_index].h -= lts_.coeff(event_index);
        kkt_matrix.aux[impulse_index].Qtt += Qts_.coeff(event_index, event_index);
        if (event_index<num_events-1) {
          kkt_matrix.aux[impulse_index].Qtt_prev += Qts_.coeff(event_index, event_index+1);
        }
        ++impulse_index;
      }
      else {
        kkt_residual.lift[lift_index].h -= lts_.coeff(event_index);
        kkt_matrix.lift[lift_index].Qtt += Qts_.coeff(event_index, event_index);
        if (event_index<num_events-1) {
          kkt_matrix.lift[lift_index].Qtt_prev += Qts_.coeff(event_index, event_index+1);
        }
        ++lift_index;
      }
    }
  }
}

} // namespace robotoc