#ifndef IDOCP_IMPULSE_COST_FUNCTION_HXX_
#define IDOCP_IMPULSE_COST_FUNCTION_HXX_

#include <cassert>

namespace idocp {

inline ImpulseCostFunction::ImpulseCostFunction()
  : costs_() {
}


inline ImpulseCostFunction::~ImpulseCostFunction() {
}


inline void ImpulseCostFunction::push_back(
    const ImpulseCostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


inline void ImpulseCostFunction::clear() {
  costs_.clear();
}


inline bool ImpulseCostFunction::isEmpty() const {
  return costs_.empty();
}


inline CostFunctionData ImpulseCostFunction::createCostFunctionData(
    const Robot& robot) const {
  auto data = CostFunctionData(robot);
  return data;
}


inline double ImpulseCostFunction::l(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->l(robot, data, t, s);
  }
  return l;
}


inline void ImpulseCostFunction::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, s, kkt_residual);
    cost->lv(robot, data, t, s, kkt_residual);
    cost->lf(robot, data, t, s, kkt_residual);
    cost->ldv(robot, data, t, s, kkt_residual);
  }
}


inline void ImpulseCostFunction::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, s, kkt_matrix);
    cost->lvv(robot, data, t, s, kkt_matrix);
    cost->lff(robot, data, t, s, kkt_matrix);
    cost->ldvdv(robot, data, t, s, kkt_matrix);
  }
}

} // namespace idocp

#endif // IDOCP_IMPULSE_COST_FUNCTION_HXX_ 