#include "robotoc/hybrid/sto_cost_function.hpp"


namespace robotoc {

STOCostFunction::STOCostFunction()
  : costs_(),
    lts_(),
    Qts_() {
}


STOCostFunction::~STOCostFunction() {
}


std::shared_ptr<STOCostFunction> STOCostFunction::clone() const {
  auto cost = std::make_shared<STOCostFunction>();
  for (const auto& e : costs_) {
    cost->push_back(e->clone());
  }
  return cost;
}


void STOCostFunction::push_back(const STOCostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


void STOCostFunction::clear() {
  costs_.clear();
}


double STOCostFunction::evalCost(const TimeDiscretization& discretization) {
  if (!costs_.empty()) {
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->evalCost(discretization);
    }
    return cost;
  }
  else {
    return 0.0;
  }
}


double STOCostFunction::linearizeCost(const TimeDiscretization& discretization, 
                                      KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = discretization.N_impulse() + discretization.N_lift();
    lts_.resize(num_events);
    lts_.setZero();
    double cost = 0;
    for (auto& e : costs_) {
      cost += e->evalCost(discretization);
      e->evalCostDerivatives(discretization, lts_);
    }
    setToKKT(discretization, kkt_residual);
    return cost;
  }
  else {
    return 0.0;
  }
}


double STOCostFunction::quadratizeCost(const TimeDiscretization& discretization, 
                                       KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual) {
  if (!costs_.empty()) {
    const int num_events = discretization.N_impulse() + discretization.N_lift();
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
      cost += e->evalCost(discretization);
      e->evalCostDerivatives(discretization, lts_);
      e->evalCostHessian(discretization, Qts_);
    }
    setToKKT(discretization, kkt_matrix, kkt_residual);
    return cost;
  }
  else {
    return 0.0;
  }
}


void STOCostFunction::setToKKT(const TimeDiscretization& discretization, 
                               KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (discretization.eventType(0) == DiscreteEventType::Impulse) {
      kkt_residual.aux[0].h -= lts_.coeff(0);
      ++impulse_index;
    }
    else {
      kkt_residual.lift[0].h -= lts_.coeff(0);
      ++lift_index;
    }
    for (int event_index=1; event_index<num_events; ++event_index) {
      assert(event_index == impulse_index+lift_index);
      if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
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


void STOCostFunction::setToKKT(const TimeDiscretization& discretization, 
                               KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  const int num_events = discretization.N_impulse() + discretization.N_lift();
  if (num_events > 0) {
    int impulse_index = 0;
    int lift_index = 0;
    if (discretization.eventType(0) == DiscreteEventType::Impulse) {
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
      if (discretization.eventType(event_index) == DiscreteEventType::Impulse) {
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