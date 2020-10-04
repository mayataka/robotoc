#ifndef IDOCP_COST_FUNCTION_HXX_
#define IDOCP_COST_FUNCTION_HXX_

#include <assert.h>

namespace idocp {

inline CostFunction::CostFunction()
  : costs_() {
}


inline CostFunction::~CostFunction() {
}


inline void CostFunction::push_back(
    const std::shared_ptr<CostFunctionComponentBase>& cost) {
  costs_.push_back(cost);
}


inline void CostFunction::clear() {
  costs_.clear();
}


inline bool CostFunction::isEmpty() const {
  return costs_.empty();
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
    KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, dtau, s, kkt_residual);
    cost->lv(robot, data, t, dtau, s, kkt_residual);
    cost->la(robot, data, t, dtau, s, kkt_residual);
    // if (robot.has_active_contacts() > 0) {
    //   cost->lf(robot, data, t, dtau, s, kkt_residual);
    // }
  }
}


inline void CostFunction::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const double dtau, const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, dtau, s, kkt_matrix);
    cost->lvv(robot, data, t, dtau, s, kkt_matrix);
    cost->laa(robot, data, t, dtau, s, kkt_matrix);
    // if (robot.has_active_contacts() > 0) {
    //   cost->lff(robot, data, t, dtau, s, kkt_matrix);
    // }
  }
}


inline void CostFunction::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, KKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->phiq(robot, data, t, s, kkt_residual);
    cost->phiv(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->phiqq(robot, data, t, s, kkt_matrix);
    cost->phivv(robot, data, t, s, kkt_matrix);
  }
}


inline void CostFunction::lq(Robot& robot, CostFunctionData& data, 
                             const double t, const double dtau, 
                             const SplitSolution& s, 
                             KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lq(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::lv(Robot& robot, CostFunctionData& data, 
                             const double t, const double dtau, 
                             const SplitSolution& s, 
                             KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lv(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::la(Robot& robot, CostFunctionData& data, 
                             const double t, const double dtau, 
                             const SplitSolution& s, 
                             KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->la(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::lf(Robot& robot, CostFunctionData& data, 
                             const double t, const double dtau, 
                             const SplitSolution& s, 
                             KKTResidual& kkt_residual) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lf(robot, data, t, dtau, s, kkt_residual);
  }
}


inline void CostFunction::lqq(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s, 
                              KKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lqq(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::lvv(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s, 
                              KKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lvv(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::laa(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s, 
                              KKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->laa(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::lff(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const SplitSolution& s, 
                              KKTMatrix& kkt_matrix) const {
  assert(dtau > 0);
  for (const auto cost : costs_) {
    cost->lff(robot, data, t, dtau, s, kkt_matrix);
  }
}


inline void CostFunction::phiq(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               KKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->phiq(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::phiv(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               KKTResidual& kkt_residual) const {
  for (const auto cost : costs_) {
    cost->phiv(robot, data, t, s, kkt_residual);
  }
}


inline void CostFunction::phiqq(Robot& robot, CostFunctionData& data, 
                                const double t, const SplitSolution& s, 
                                KKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->phiqq(robot, data, t, s, kkt_matrix);
  }
}


inline void CostFunction::phivv(Robot& robot, CostFunctionData& data, 
                                const double t, const SplitSolution& s, 
                                KKTMatrix& kkt_matrix) const {
  for (const auto cost : costs_) {
    cost->phivv(robot, data, t, s, kkt_matrix);
  }
}


inline void CostFunction::lu(Robot& robot, CostFunctionData& data, 
                             const double t, const double dtau, 
                             const Eigen::VectorXd& u, 
                             Eigen::VectorXd& lu) const {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(lu.size() == robot.dimv());
  for (const auto cost : costs_) {
    cost->lu(robot, data, t, dtau, u, lu);
  }
}


inline void CostFunction::luu(Robot& robot, CostFunctionData& data, 
                              const double t, const double dtau, 
                              const Eigen::VectorXd& u, 
                              Eigen::MatrixXd& Quu) const {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  assert(Quu.rows() == robot.dimv());
  assert(Quu.cols() == robot.dimv());
  for (const auto cost : costs_) {
    cost->luu(robot, data, t, dtau, u, Quu);
  }
}

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_HXX_