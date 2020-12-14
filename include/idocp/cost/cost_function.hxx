#ifndef IDOCP_COST_FUNCTION_HXX_
#define IDOCP_COST_FUNCTION_HXX_

#include <cassert>

namespace idocp {

inline CostFunction::CostFunction()
  : costs_(),
    impulse_cost_function_(std::make_shared<ImpulseCostFunction>()) {
}


inline CostFunction::~CostFunction() {
}


inline void CostFunction::push_back(const CostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


inline void CostFunction::push_back(
    const ImpulseCostFunctionComponentBasePtr& cost) {
  impulse_cost_function_->push_back(cost);
}


inline std::shared_ptr<ImpulseCostFunction> 
CostFunction::getImpulseCostFunction() {
  return impulse_cost_function_;
}


inline void CostFunction::clear() {
  costs_.clear();
  impulse_cost_function_->clear();
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


inline double CostFunction::l(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s) const {
  assert(dtau > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->l(robot, data, t, dtau, s);
  }
  return l;
}


inline double CostFunction::phi(Robot& robot, CostFunctionData& data, 
                                const double t, const SplitSolution& s) const {
  double phi = 0;
  for (const auto cost : costs_) {
    phi += cost->phi(robot, data, t, s);
  }
  return phi;
}


inline void CostFunction::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const double dtau, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, dtau, s, kkt_residual);
    cost->lv(robot, data, t, dtau, s, kkt_residual);
    cost->la(robot, data, t, dtau, s, kkt_residual);
    if (s.dimf() > 0) {
      cost->lf(robot, data, t, dtau, s, kkt_residual);
    }
    cost->lu(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const double dtau, const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, dtau, s, kkt_matrix);
    cost->lvv(robot, data, t, dtau, s, kkt_matrix);
    cost->laa(robot, data, t, dtau, s, kkt_matrix);
    if (s.dimf() > 0) {
      cost->lff(robot, data, t, dtau, s, kkt_matrix);
    }
    cost->luu(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->phiq(robot, data, t, s, kkt_residual);
    cost->phiv(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->phiqq(robot, data, t, s, kkt_matrix);
    cost->phivv(robot, data, t, s, kkt_matrix);
  }
}

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_HXX_