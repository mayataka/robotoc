#include "robotoc/cost/cost_function.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

CostFunction::CostFunction(const double discount_factor)
  : costs_(),
    discount_factor_(discount_factor),
    discounted_cost_(true) {
  try {
    if (discount_factor <= 0.0) {
      throw std::out_of_range("invalid argument: discount_factor must be positive!");
    }
    if (discount_factor >= 1.0) {
      throw std::out_of_range("invalid argument: discount_factor must be smaller than 1.0!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


CostFunction::CostFunction()
  : costs_(),
    discount_factor_(1.0),
    discounted_cost_(false) {
}


CostFunction::~CostFunction() {
}


void CostFunction::setDiscountFactor(const double discount_factor) {
  if (discount_factor > 0 && discount_factor < 1.0) {
    discount_factor_ = discount_factor;
    discounted_cost_ = true;
  }
  else {
    discount_factor_ = 1.0;
    discounted_cost_ = false;
  }
}


double CostFunction::discountFactor() const {
  if (discounted_cost_) {
    return discount_factor_;
  }
  else {
    return 1.0;
  }
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
    kkt_residual.lu.array() *= f;
    kkt_residual.la.array() *= f;
    kkt_residual.lf().array() *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
    kkt_residual.lu.array() *= f;
    kkt_residual.la.array() *= f;
    kkt_residual.lf().array() *= f;
    kkt_matrix.Qxx.array() *= f;
    kkt_matrix.Qxu.array() *= f;
    kkt_matrix.Quu.array() *= f;
    kkt_matrix.Qaa.array() *= f;
    kkt_matrix.Qff().array() *= f;
    kkt_matrix.Qqf().array() *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
    kkt_matrix.Qxx.array() *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
    kkt_residual.ldv.array() *= f;
    kkt_residual.lf().array() *= f;
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
  if (discounted_cost_) {
    const double f = pow(discount_factor_, grid_info.time_stage);
    l *= f;
    kkt_residual.lx.array() *= f;
    kkt_residual.ldv.array() *= f;
    kkt_residual.lf().array() *= f;
    kkt_matrix.Qxx.array() *= f;
    kkt_matrix.Qdvdv.array() *= f;
    kkt_matrix.Qff().array() *= f;
    kkt_matrix.Qqf().array() *= f;
  }
  return l;
}

} // namespace robotoc