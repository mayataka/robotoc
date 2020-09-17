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

///
/// @class CostFunction
/// @brief Cost function.
///
class CostFunction {
public:

  ///
  /// @brief Default constructor. 
  ///
  CostFunction();

  ///
  /// @brief Destructor. 
  ///
  ~CostFunction();

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunction(const CostFunction&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunction& operator=(const CostFunction&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunction(CostFunction&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunction& operator=(CostFunction&&) noexcept = default;

  ///
  /// @brief Append a cost function component to the cost function.
  /// @param[in] cost Cost function component appended to the cost.
  ///
  void push_back(const std::shared_ptr<CostFunctionComponentBase>& cost);

  ///
  /// @brief Clear cost function by removing all components.
  ///
  void clear();

  ///
  /// @brief Check whether the cost function is empty or not.
  /// @return true if the cost function is empty. false if not.
  ///
  bool isEmpty() const;

  ///
  /// @brief Return true if the cost function component requres kinematics of 
  /// robot model. Return false if not.
  /// @return true if the cost function component requres kinematics of 
  /// robot model. false if not.
  ///
  bool useKinematics() const;

  ///
  /// @brief Creates CostFunctionData. 
  /// @return Cost function data.
  ///
  CostFunctionData createCostFunctionData(const Robot& robot) const;

  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const;

  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const;

  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   KKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               KKTMatrix& kkt_matrix) const;

  void computeTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      KKTResidual& kkt_residual) const;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  KKTMatrix& kkt_matrix) const;

  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  void phiq(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  void phiv(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const;

  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const double dtau, const SplitSolution& s, 
             KKTMatrix& kkt_matrix) const;

  void lu(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const;

  void luu(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const;

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp

#include "idocp/cost/cost_function.hxx"

#endif // IDOCP_COST_FUNCTION_HPP_