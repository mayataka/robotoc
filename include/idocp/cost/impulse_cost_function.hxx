#ifndef IDOCP_IMPULSE_COST_FUNCTION_HXX_
#define IDOCP_IMPULSE_COST_FUNCTION_HXX_

#include <assert.h>

namespace idocp {

inline ImpulseCostFunction::ImpulseCostFunction()
  : costs_() {
}


inline ImpulseCostFunction::~ImpulseCostFunction() {
}


inline void ImpulseCostFunction::push_back(
    const std::shared_ptr<ImpulseCostFunctionComponentBase>& cost) {
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
    const ImpulseSplitSolution& s, ImpulseKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, s, kkt_residual);
    cost->lv(robot, data, t, s, kkt_residual);
    // if (robot.has_active_contacts() > 0) {
    //   cost->lf(robot, data, t, s, kkt_residual);
    // }
  }
}


inline void ImpulseCostFunction::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, s, kkt_matrix);
    cost->lvv(robot, data, t, s, kkt_matrix);
    // if (robot.has_active_contacts() > 0) {
    //   cost->lff(robot, data, t, s, kkt_matrix);
    // }
  }
}


inline void ImpulseCostFunction::lq(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const ImpulseSplitSolution& s, 
                                    ImpulseKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, s, kkt_residual);
  }
}


inline void ImpulseCostFunction::lv(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const ImpulseSplitSolution& s, 
                                    ImpulseKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->lv(robot, data, t, s, kkt_residual);
  }
}


inline void ImpulseCostFunction::lf(Robot& robot, CostFunctionData& data, 
                                    const double t, 
                                    const ImpulseSplitSolution& s, 
                                    ImpulseKKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->lf(robot, data, t, s, kkt_residual);
  }
}


inline void ImpulseCostFunction::lqq(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, s, kkt_matrix);
  }
}


inline void ImpulseCostFunction::lvv(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->lvv(robot, data, t, s, kkt_matrix);
  }
}


inline void ImpulseCostFunction::lff(Robot& robot,  CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseKKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->lff(robot, data, t, s, kkt_matrix);
  }
}


inline void ImpulseCostFunction::ldv(Robot& robot, CostFunctionData& data, 
                                     const double t, const Eigen::VectorXd& dv, 
                                     Eigen::VectorXd& ldv) const {
  assert(dv.size() == robot.dimv());
  assert(ldv.size() == robot.dimv());
  for (const auto cost : costs_) {
    cost->ldv(robot, data, t, dv, ldv);
  }
}


inline void ImpulseCostFunction::ldvdv(Robot& robot, CostFunctionData& data, 
                                       const double t, 
                                       const Eigen::VectorXd& dv, 
                                       Eigen::MatrixXd& Qdvdv) const {
  assert(dv.size() == robot.dimv());
  assert(Qdvdv.rows() == robot.dimv());
  assert(Qdvdv.cols() == robot.dimv());
  for (const auto cost : costs_) {
    cost->ldvdv(robot, data, t, dv, Qdvdv);
  }
}

} // namespace idocp

#endif // IDOCP_IMPULSE_COST_FUNCTION_HXX_ 