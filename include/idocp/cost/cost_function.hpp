#ifndef IDOCP_COST_FUNCTION_HPP_
#define IDOCP_COST_FUNCTION_HPP_

#include <vector>
#include <memory>
#include <assert.h>

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class CostFunction {
public:
  CostFunction() 
    : costs_() {
  }

  ~CostFunction() {}

  // Use default copy constructor.
  CostFunction(const CostFunction&) = default;

  // Use default copy coperator.
  CostFunction& operator=(const CostFunction&) = default;

  // Use default move constructor.
  CostFunction(CostFunction&&) noexcept = default;

  // Use default move assign coperator.
  CostFunction& operator=(CostFunction&&) noexcept = default;

  void push_back(const std::shared_ptr<CostFunctionComponentBase>& cost) {
    costs_.push_back(cost);
  }

  void clear() {
    costs_.clear();
  }

  bool isEmpty() {
    return costs_.empty();
  }

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const {
    assert(dtau > 0);
    double l = 0;
    for (const auto cost : costs_) {
      l += cost->l(robot, data, t, dtau, s);
    }
    return l;
  }

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const {
    double phi = 0;
    for (const auto cost : costs_) {
      phi += cost->phi(robot, data, t, s);
    }
    return phi;
  }

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lq(robot, data, t, dtau, s, kkt_residual);
    }
  }

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lv(robot, data, t, dtau, s, kkt_residual);
    }
  }

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->la(robot, data, t, dtau, s, kkt_residual);
    }
  }

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lf(robot, data, t, dtau, s, kkt_residual);
    }
  }

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lu(robot, data, t, dtau, s, kkt_residual);
    }
  }

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lqq(robot, data, t, dtau, s, kkt_matrix);
    }
  }

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lvv(robot, data, t, dtau, s, kkt_matrix);
    }
  }

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->laa(robot, data, t, dtau, s, kkt_matrix);
    }
  }

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->lff(robot, data, t, dtau, s, kkt_matrix);
    }
  }

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const {
    assert(dtau > 0);
    for (const auto cost : costs_) {
      cost->luu(robot, data, t, dtau, s, kkt_matrix);
    }
  }

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const {
    for (const auto cost : costs_) {
      cost->phiq(robot, data, t, s, kkt_residual);
    }
  }

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const {
    for (const auto cost : costs_) {
      cost->phiv(robot, data, t, s, kkt_residual);
    }
  }

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const {
    for (const auto cost : costs_) {
      cost->phiqq(robot, data, t, s, kkt_matrix);
    }
  }

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const {
    for (const auto cost : costs_) {
      cost->phivv(robot, data, t, s, kkt_matrix);
    }
  }

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_HPP_