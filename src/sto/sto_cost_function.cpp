#include "robotoc/sto/sto_cost_function.hpp"


namespace robotoc {

STOCostFunction::STOCostFunction()
  : costs_() {
}


bool STOCostFunction::exist(const std::string& name) const {
  return (cost_names_.find(name) != cost_names_.end());
}


void STOCostFunction::add(const std::string& name, 
                          const STOCostFunctionComponentBasePtr& cost) {
  if (exist(name)) {
    throw std::runtime_error("[STOCostFunction] invalid argument: cost component '" + name + "' already exists!");
  }
  cost_names_.emplace(name, costs_.size());
  costs_.push_back(cost);
}


void STOCostFunction::erase(const std::string& name) {
  if (!exist(name)) {
    throw std::runtime_error("[STOCostFunction] invalid argument: cost component '" + name + "' does not exist!");
  }
  const int index = cost_names_.at(name);
  cost_names_.erase(name);
  costs_.erase(costs_.begin()+index);
}


std::shared_ptr<STOCostFunctionComponentBase> 
STOCostFunction::get(const std::string& name) const {
  if (!exist(name)) {
    throw std::runtime_error("[STOCostFunction] invalid argument: cost component '" + name + "' does not exist!");
  }
  const int index = cost_names_.at(name);
  return costs_.at(index);
}


void STOCostFunction::clear() {
  costs_.clear();
  cost_names_.clear();
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