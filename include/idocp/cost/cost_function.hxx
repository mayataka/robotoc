#ifndef IDOCP_COST_FUNCTION_HXX_
#define IDOCP_COST_FUNCTION_HXX_

#include <cassert>

namespace idocp {

inline CostFunction::CostFunction()
  : costs_() {
}


inline CostFunction::~CostFunction() {
}


inline void CostFunction::push_back(const CostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


inline void CostFunction::clear() {
  costs_.clear();
}


inline bool CostFunction::useKinematics() const {
  for (const auto cost : costs_) {
    if (cost->useKinematics()) {
      return true;
    }
  }
  return false;
}


inline CostFunctionData CostFunction::createCostFunctionData(
    const Robot& robot) const {
  auto data = CostFunctionData(robot);
  return data;
}


inline double CostFunction::computeStageCost(
    Robot& robot,  CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s) const {
  assert(dtau >= 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeStageCost(robot, data, t, dtau, s);
  }
  return l;
}


inline double CostFunction::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeTerminalCost(robot, data, t, s);
  }
  return l;
}


inline double CostFunction::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeImpulseCost(robot, data, t, s);
  }
  return l;
}


inline void CostFunction::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dtau >= 0);
  for (const auto cost : costs_) {
    cost->computeStageCostDerivatives(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  assert(dtau >= 0);
  for (const auto cost : costs_) {
    cost->computeStageCostHessian(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->computeTerminalCostHessian(robot, data, t, s, kkt_matrix);
  }
}


inline void CostFunction::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->computeImpulseCostHessian(robot, data, t, s, kkt_matrix);
  }
}

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_HXX_