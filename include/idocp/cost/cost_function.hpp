#ifndef IDOCP_COST_FUNCTION_HPP_
#define IDOCP_COST_FUNCTION_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class CostFunction {
public:
  CostFunction();

  ~CostFunction();

  // Use default copy constructor.
  CostFunction(const CostFunction&) = default;

  // Use default copy coperator.
  CostFunction& operator=(const CostFunction&) = default;

  // Use default move constructor.
  CostFunction(CostFunction&&) noexcept = default;

  // Use default move assign coperator.
  CostFunction& operator=(CostFunction&&) noexcept = default;

  void push_back(const std::shared_ptr<CostFunctionComponentBase>& cost);

  void clear();

  bool isEmpty() const;

  double l(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const;

  double phi(const Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const;

  void computeStageCostDerivatives(const Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   KKTResidual& kkt_residual) const;

  void computeStageCostHessian(const Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               KKTMatrix& kkt_matrix) const;

  void computeTerminalCostDerivatives(const Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      KKTResidual& kkt_residual) const;

  void computeTerminalCostHessian(const Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  KKTMatrix& kkt_matrix) const;

  void lq(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lv(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void la(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lf(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lqq(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void lvv(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void laa(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void lff(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void phiq(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  void phiv(const Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const;

  void phivv(const Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const;

  void lu(const Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const;

  void luu(const Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const;

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp

#include "idocp/cost/cost_function.hxx"

#endif // IDOCP_COST_FUNCTION_HPP_