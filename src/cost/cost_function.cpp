#include "robotoc/cost/cost_function.hpp"

#include <cassert>


namespace robotoc {

CostFunction::CostFunction()
  : costs_() {
}


CostFunction::~CostFunction() {
}


void CostFunction::push_back(const CostFunctionComponentBasePtr& cost) {
  costs_.push_back(cost);
}


void CostFunction::clear() {
  costs_.clear();
}


bool CostFunction::useKinematics() const {
  for (const auto cost : costs_) {
    if (cost->useKinematics()) {
      return true;
    }
  }
  return false;
}


CostFunctionData CostFunction::createCostFunctionData(const Robot& robot) const {
  auto data = CostFunctionData(robot);
  return data;
}


double CostFunction::evalStageCost(Robot& robot, 
                                   const ContactStatus& contact_status, 
                                   CostFunctionData& data, 
                                   const GridInfo& grid_info,
                                   const SplitSolution& s) const {
  assert(grid_info.dt > 0);
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalStageCost(robot, contact_status, data, grid_info, s);
  }
  return l;
}


double CostFunction::linearizeStageCost(Robot& robot, 
                                        const ContactStatus& contact_status, 
                                        CostFunctionData& data, 
                                        const GridInfo& grid_info, 
                                        const SplitSolution& s, 
                                        SplitKKTResidual& kkt_residual) const {
  assert(grid_info.dt > 0);
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalStageCost(robot, contact_status, data, grid_info, s);
    e->evalStageCostDerivatives(robot, contact_status, data, grid_info, s,
                                kkt_residual);
  }
  return l;
}


double CostFunction::quadratizeStageCost(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         CostFunctionData& data, 
                                         const GridInfo& grid_info, 
                                         const SplitSolution& s, 
                                         SplitKKTResidual& kkt_residual, 
                                         SplitKKTMatrix& kkt_matrix) const {
  assert(grid_info.dt > 0);
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalStageCost(robot, contact_status, data, grid_info, s);
    e->evalStageCostDerivatives(robot, contact_status, data, grid_info, s,
                                kkt_residual);
    e->evalStageCostHessian(robot, contact_status, data, grid_info, s,
                            kkt_matrix);
  }
  return l;
}


double CostFunction::evalTerminalCost(Robot& robot, CostFunctionData& data, 
                                      const GridInfo& grid_info, 
                                      const SplitSolution& s) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalTerminalCost(robot, data, grid_info, s);
  }
  return l;
}


double CostFunction::linearizeTerminalCost(Robot& robot, CostFunctionData& data, 
                                           const GridInfo& grid_info, 
                                           const SplitSolution& s, 
                                           SplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalTerminalCost(robot, data, grid_info, s);
    e->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_residual);
  }
  return l;
}


double CostFunction::quadratizeTerminalCost(Robot& robot, 
                                            CostFunctionData& data, 
                                            const GridInfo& grid_info, 
                                            const SplitSolution& s, 
                                            SplitKKTResidual& kkt_residual, 
                                            SplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalTerminalCost(robot, data, grid_info, s);
    e->evalTerminalCostDerivatives(robot, data, grid_info, s, kkt_residual);
    e->evalTerminalCostHessian(robot, data, grid_info, s, kkt_matrix);
  }
  return l;
}


double CostFunction::evalImpulseCost(Robot& robot, 
                                     const ImpulseStatus& impulse_status, 
                                     CostFunctionData& data, 
                                     const GridInfo& grid_info, 
                                     const ImpulseSplitSolution& s) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalImpulseCost(robot, impulse_status, data, grid_info, s);
  }
  return l;
}


double CostFunction::linearizeImpulseCost(Robot& robot, 
                                          const ImpulseStatus& impulse_status, 
                                          CostFunctionData& data, 
                                          const GridInfo& grid_info, 
                                          const ImpulseSplitSolution& s, 
                                          ImpulseSplitKKTResidual& kkt_residual) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalImpulseCost(robot, impulse_status, data, grid_info, s);
    e->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, 
                                  kkt_residual);
  }
  return l;
}


double CostFunction::quadratizeImpulseCost(Robot& robot, 
                                           const ImpulseStatus& impulse_status, 
                                           CostFunctionData& data, 
                                           const GridInfo& grid_info, 
                                           const ImpulseSplitSolution& s, 
                                           ImpulseSplitKKTResidual& kkt_residual, 
                                           ImpulseSplitKKTMatrix& kkt_matrix) const {
  double l = 0;
  for (const auto e : costs_) {
    l += e->evalImpulseCost(robot, impulse_status, data, grid_info, s);
    e->evalImpulseCostDerivatives(robot, impulse_status, data, grid_info, s, 
                                  kkt_residual);
    e->evalImpulseCostHessian(robot, impulse_status, data, grid_info, s, 
                              kkt_matrix);
  }
  return l;
}

} // namespace robotoc