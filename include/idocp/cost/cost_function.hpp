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
/// @brief Stack of the cost function. Composed by cost function components 
/// that inherits CostFunctionComponentBase.
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
  /// @param[in] cost shared pointer to the cost function component appended 
  /// to the cost.
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
  /// @brief Check if the cost function component requres kinematics of robot 
  /// model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  bool useKinematics() const;

  ///
  /// @brief Creates CostFunctionData according to robot model and cost 
  /// function components. 
  /// @param[in] robot robot model.
  /// @return Cost function data.
  ///
  CostFunctionData createCostFunctionData(const Robot& robot) const;

  ///
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  double l(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s) const;

  ///
  /// @brief Computes and returns terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  double phi(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration, velocity, acceleration, and contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dtau, 
                               const SplitSolution& s, 
                               KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the partial derivatives of the terminal cost with respect
  /// to the configuration and velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void computeTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the configuration and velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void lq(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void lv(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the acceleration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void la(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void lf(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const SplitSolution& s, 
          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void lqq(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void lvv(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the acceleration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void laa(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void lff(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const SplitSolution& s, 
           KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the partial derivatives of the terminal cost with respect
  /// to the configuration.
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void phiq(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the partial derivatives of the terminal cost with respect
  /// to the velocity.
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual The KKT residual. The partial derivatives are 
  /// added to this data.
  ///
  void phiv(Robot& robot, CostFunctionData& data, const double t, 
            const SplitSolution& s, KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void phiqq(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  void phivv(Robot& robot, CostFunctionData& data, const double t, 
             const SplitSolution& s, KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Computes the partial derivatives of the stage cost with respect
  /// to the control input torques. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] u Control input torques. Size must be Robot::dimv().
  /// @param[out] lu The KKT residual with respect to u. Size must be 
  /// Robot::dimv().
  ///
  void lu(Robot& robot, CostFunctionData& data, const double t, 
          const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu) const;

  ///
  /// @brief Computes the Hessian of the stage cost with respect
  /// to the control input torques. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] u Control input torques. Size must be Robot::dimv().
  /// @param[out] Quu The Hessian of the KKT residual with respect to u.  
  /// Size must be Robot::dimv() x Robot::dimv().
  ///
  void luu(Robot& robot, CostFunctionData& data, const double t, 
           const double dtau, const Eigen::VectorXd& u, 
           Eigen::MatrixXd& Quu) const;

private:
  std::vector<std::shared_ptr<CostFunctionComponentBase>> costs_;

};

} // namespace idocp

#include "idocp/cost/cost_function.hxx"

#endif // IDOCP_COST_FUNCTION_HPP_