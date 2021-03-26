#ifndef IDOCP_CONFIGURATION_SPACE_COST_HPP_
#define IDOCP_CONFIGURATION_SPACE_COST_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

///
/// @class ConfigurationSpaceCost
/// @brief Configuration space quadratic cost. 
///
class ConfigurationSpaceCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  ConfigurationSpaceCost(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ConfigurationSpaceCost();

  ///
  /// @brief Destructor. 
  ///
  ~ConfigurationSpaceCost();

  ///
  /// @brief Default copy constructor. 
  ///
  ConfigurationSpaceCost(const ConfigurationSpaceCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConfigurationSpaceCost& operator=(const ConfigurationSpaceCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConfigurationSpaceCost(ConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConfigurationSpaceCost& operator=(ConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Sets the reference configuration q. 
  /// @param[in] q_ref Reference configuration q. Size must be Robot::dimq().
  ///
  void set_q_ref(const Eigen::VectorXd& q_ref);

  ///
  /// @brief Sets the reference velocity v. 
  /// @param[in] q_ref Reference velocity v. Size must be Robot::dimv().
  ///
  void set_v_ref(const Eigen::VectorXd& v_ref);

  ///
  /// @brief Sets the reference control input torques u. 
  /// @param[in] q_ref Reference control input torques u. Size must be 
  /// Robot::dimu().
  ///
  void set_u_ref(const Eigen::VectorXd& u_ref);

  ///
  /// @brief Sets the weight vector on the configuration q. 
  /// @param[in] q_weight Weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  ///
  void set_q_weight(const Eigen::VectorXd& q_weight);

  ///
  /// @brief Sets the weight on the velocity v. 
  /// @param[in] v_weight Weight vector on the velocity v. 
  /// Size must be Robot::dimv().
  ///
  void set_v_weight(const Eigen::VectorXd& v_weight);

  ///
  /// @brief Sets the weight on the acceleration a. 
  /// @param[in] a_weight Weight vector on the acceleration a. 
  /// Size must be Robot::dimv().
  ///
  void set_a_weight(const Eigen::VectorXd& a_weight);

  ///
  /// @brief Sets the weight on the control input torques u. 
  /// @param[in] u_weight Weight vector on the control input torques u. 
  /// Size must be Robot::dimu().
  ///
  void set_u_weight(const Eigen::VectorXd& u_weight);

  ///
  /// @brief Sets the terminal weight vector on the configuration q. 
  /// @param[in] qf_weight Terminal weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  ///
  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  ///
  /// @brief Sets the terminal weight vector on the velocity v. 
  /// @param[in] vf_weight Terminal weight vector on the velocity v. 
  /// Size must be Robot::dimv().
  ///
  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  ///
  /// @brief Sets the weight vector on the configuration q at impulse. 
  /// @param[in] qi_weight Weight vector on the configuration q at impulse. 
  /// Size must be Robot::dimv().
  ///
  void set_qi_weight(const Eigen::VectorXd& qi_weight);

  ///
  /// @brief Sets the weight vector on the velocity v at impulse. 
  /// @param[in] vi_weight Weight vector on the velocity v at impulse. 
  /// Size must be Robot::dimv().
  ///
  void set_vi_weight(const Eigen::VectorXd& vi_weight);

  ///
  /// @brief Sets the weight vector on the impulse change in the velocity dv. 
  /// @param[in] dvi_weight Weight vector on the impulse change in the velocity. 
  /// Size must be Robot::dimv().
  ///
  void set_dvi_weight(const Eigen::VectorXd& dvi_weight);

  bool useKinematics() const override;

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, const SplitSolution& s) const;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, const SplitSolution& s) const;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, const ImpulseSplitSolution& s) const;

  void computeStageCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dt, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const;

  void computeTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const;

  void computeImpulseCostDerivatives(Robot& robot, CostFunctionData& data, 
                                     const double t, 
                                     const ImpulseSplitSolution& s, 
                                     ImpulseSplitKKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const;

  void computeImpulseCostHessian(Robot& robot, CostFunctionData& data, 
                                 const double t, const ImpulseSplitSolution& s, 
                                 ImpulseSplitKKTMatrix& kkt_matrix) const;

private:
  int dimq_, dimv_, dimu_;
  Eigen::VectorXd q_ref_, v_ref_, u_ref_, q_weight_, v_weight_, a_weight_, 
                  u_weight_, qf_weight_, vf_weight_, qi_weight_, vi_weight_, 
                  dvi_weight_;
};

} // namespace idocp


#endif // IDOCP_CONFIGURATION_SPACE_COST_HPP_ 