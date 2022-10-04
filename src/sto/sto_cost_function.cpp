#include "robotoc/sto/sto_cost_function.hpp"


namespace robotoc {

STOCostFunction::STOCostFunction()
  : costs_() {
}


void STOCostFunction::push_back(const STOCostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


void STOCostFunction::clear() {
  costs_.clear();
}


double STOCostFunction::evalCost(
    const TimeDiscretization& time_discretization) const {
  if (costs_.empty()) return 0.0;

  double cost = 0;
  for (auto& e : costs_) {
    cost += e->evalCost(time_discretization);
  }
  return cost;
}


double STOCostFunction::linearizeCost(const TimeDiscretization& time_discretization, 
                                      Eigen::VectorXd& lt) const {
  if (costs_.empty()) return 0.0;

  double cost = 0;
  for (auto& e : costs_) {
    cost += e->evalCost(time_discretization);
    e->evalCostDerivatives(time_discretization, lt);
  }
  return cost;
}


double STOCostFunction::quadratizeCost(const TimeDiscretization& time_discretization, 
                                       Eigen::VectorXd& lt, 
                                       Eigen::MatrixXd& Qtt) const {
  if (costs_.empty()) return 0.0;

  double cost = 0;
  for (auto& e : costs_) {
    cost += e->evalCost(time_discretization);
    e->evalCostDerivatives(time_discretization, lt);
    e->evalCostHessian(time_discretization, Qtt);
  }
  return cost;
}

} // namespace robotoc