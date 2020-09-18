#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

///
/// @class CostFunctionComponentBase
/// @brief Base class of components of cost function.
///
class CostFunctionComponentBase {
public:

  ///
  /// @brief Default constructor. 
  ///
  CostFunctionComponentBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~CostFunctionComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Check if the cost function component requres kinematics of robot 
  /// model.
  /// @return true if the cost function component requres kinematics of 
  /// Robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  ///
  /// @brief Computes and returns stage cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @return Stage cost.
  ///
  virtual double l(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s) const = 0;

  ///
  /// @brief Computes and returns terminal cost. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @return Terminal cost.
  ///
  virtual double phi(Robot& robot, CostFunctionData& data, const double t, 
                     const SplitSolution& s) const = 0;

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
  virtual void lq(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

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
  virtual void lv(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

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
  virtual void la(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

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
  virtual void lf(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const SplitSolution& s, 
                  KKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void lqq(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void lvv(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the acceleration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void laa(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the stage cost with respect
  /// to the contact forces. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] dtau Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void lff(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const SplitSolution& s, 
                   KKTMatrix& kkt_matrix) const = 0;

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
  virtual void phiq(Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

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
  virtual void phiv(Robot& robot, CostFunctionData& data,  
                    const double t, const SplitSolution& s, 
                    KKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void phiqq(Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

  ///
  /// @brief Computes the Hessians of the terminal cost with respect
  /// to the velocity. 
  /// @param[in] robot Robot model.
  /// @param[in] data Cost function data.
  /// @param[in] t Time.
  /// @param[in] s Split solution.
  /// @param[out] kkT_matrix The KKT matrix. The Hessians are added to this 
  /// data.
  ///
  virtual void phivv(Robot& robot, CostFunctionData& data,  
                     const double t, const SplitSolution& s,
                     KKTMatrix& kkt_matrix) const = 0;

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
  virtual void lu(Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lu) const = 0;

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
  virtual void luu(Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& Quu) const  = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_