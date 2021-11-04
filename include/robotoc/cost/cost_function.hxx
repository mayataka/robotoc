#ifndef ROBOTOC_COST_FUNCTION_HXX_
#define ROBOTOC_COST_FUNCTION_HXX_

#include <cassert>

namespace robotoc {

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


inline double CostFunction::evalStageCost(Robot& robot, 
                                          const ContactStatus& contact_status, 
                                          CostFunctionData& data, 
                                          const double t, const double dt, 
                                          const SplitSolution& s) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalStageCost(robot, contact_status, data, t, dt, s);
  }
  return l;
}


inline double CostFunction::linearizeStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalStageCost(robot, contact_status, data, t, dt, s);
    cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, 
                                   kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const double t, const double dt, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual, SplitKKTMatrix& kkt_matrix) const {
  assert(dt > 0);
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalStageCost(robot, contact_status, data, t, dt, s);
    cost->evalStageCostDerivatives(robot, contact_status, data, t, dt, s, 
                                   kkt_residual);
    cost->evalStageCostHessian(robot, contact_status, data, t, dt, s, 
                               kkt_matrix);
  }
  return l;
}


inline double CostFunction::evalTerminalCost(Robot& robot, 
                                             CostFunctionData& data, 
                                             const double t, 
                                             const SplitSolution& s) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalTerminalCost(robot, data, t, s);
  }
  return l;
}


inline double CostFunction::linearizeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalTerminalCost(robot, data, t, s);
    cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual, 
    SplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalTerminalCost(robot, data, t, s);
    cost->evalTerminalCostDerivatives(robot, data, t, s, kkt_residual);
    cost->evalTerminalCostHessian(robot, data, t, s, kkt_matrix);
  }
  return l;
}


inline double CostFunction::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalImpulseCost(robot, impulse_status, data, t, s);
  }
  return l;
}


inline double CostFunction::linearizeImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalImpulseCost(robot, impulse_status, data, t, s);
    cost->evalImpulseCostDerivatives(robot, impulse_status, data, t, s, 
                                     kkt_residual);
  }
  return l;
}


inline double CostFunction::quadratizeImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const double t, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto cost : costs_) {
    l += cost->evalImpulseCost(robot, impulse_status, data, t, s);
    cost->evalImpulseCostDerivatives(robot, impulse_status, data, t, s, 
                                     kkt_residual);
    cost->evalImpulseCostHessian(robot, impulse_status, data, t, s, kkt_matrix);
  }
  return l;
}

} // namespace robotoc

#endif // ROBOTOC_COST_FUNCTION_HXX_