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
    Robot& robot,  CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeStageCost(robot, data, t, dt, s);
  }
  return l;
}


inline double CostFunction::linearizeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeStageCost(robot, data, t, dt, s);
    cost->computeStageCostDerivatives(robot, data, t, dt, s, kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual,
    SplitKKTMatrix& kkt_matrix) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeStageCost(robot, data, t, dt, s);
    cost->computeStageCostDerivatives(robot, data, t, dt, s, kkt_residual);
    cost->computeStageCostHessian(robot, data, t, dt, s, kkt_matrix);
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


inline double CostFunction::linearizeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeTerminalCost(robot, data, t, s);
    cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeTerminalCost(robot, data, t, s);
    cost->computeTerminalCostDerivatives(robot, data, t, s, kkt_residual);
    cost->computeTerminalCostHessian(robot, data, t, s, kkt_matrix);
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


inline double CostFunction::linearizeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeImpulseCost(robot, data, t, s);
    cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual,
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->computeImpulseCost(robot, data, t, s);
    cost->computeImpulseCostDerivatives(robot, data, t, s, kkt_residual);
    cost->computeImpulseCostHessian(robot, data, t, s, kkt_matrix);
  }
  return l;
}

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_HXX_