#ifndef ROBOTOC_CONFIGURATION_SPACE_COST_HPP_
#define ROBOTOC_CONFIGURATION_SPACE_COST_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class ConfigurationSpaceCost
/// @brief Configuration space cost. 
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
  /// @param[in] v_ref Reference velocity v. Size must be Robot::dimv().
  ///
  void set_v_ref(const Eigen::VectorXd& v_ref);

  ///
  /// @brief Sets the reference control input torques u. 
  /// @param[in] u_ref Reference control input torques u. Size must be 
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
  /// @brief Sets the weight vector on the configuration q at the terminal stage. 
  /// @param[in] qf_weight Weight vector on the configuration q at the terminal 
  /// stage. Size must be Robot::dimv().
  ///
  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  ///
  /// @brief Sets the weight vector on the velocity v at the terminal stage. 
  /// @param[in] vf_weight Weight vector on the velocity v at the terminal 
  /// stage. Size must be Robot::dimv().
  ///
  void set_vf_weight(const Eigen::VectorXd& vf_weight);

  ///
  /// @brief Sets the weight vector on the configuration q at impulse stages. 
  /// @param[in] qi_weight Weight vector on the configuration q at impulse  
  /// stages. Size must be Robot::dimv().
  ///
  void set_qi_weight(const Eigen::VectorXd& qi_weight);

  ///
  /// @brief Sets the weight vector on the velocity v at the impulse stages. 
  /// @param[in] vi_weight Weight vector on the velocity v at the impulse  
  /// stages. Size must be Robot::dimv().
  ///
  void set_vi_weight(const Eigen::VectorXd& vi_weight);

  ///
  /// @brief Sets the weight vector on the impulse change in the velocity dv at 
  /// the impulse stages. 
  /// @param[in] dvi_weight Weight vector on the impulse change in the velocity
  /// the impulse stages. Size must be Robot::dimv().
  ///
  void set_dvi_weight(const Eigen::VectorXd& dvi_weight);

  bool useKinematics() const override;

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const GridInfo& grid_info, 
                       const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const GridInfo& grid_info, 
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const GridInfo& grid_info, 
                            const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const GridInfo& grid_info, const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const GridInfo& grid_info, const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const GridInfo& grid_info, const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const GridInfo& grid_info, 
                         const ImpulseSplitSolution& s) const override;

  void evalImpulseCostDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                                  CostFunctionData& data, const GridInfo& grid_info, 
                                  const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const override;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const GridInfo& grid_info, 
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

private:
  int dimq_, dimv_, dimu_;
  Eigen::VectorXd q_ref_, v_ref_, u_ref_, q_weight_, v_weight_, a_weight_, 
                  u_weight_, qf_weight_, vf_weight_, qi_weight_, vi_weight_, 
                  dvi_weight_;
};

} // namespace robotoc


#endif // ROBOTOC_CONFIGURATION_SPACE_COST_HPP_ 