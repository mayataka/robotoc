#ifndef ROBOTOC_MULTI_MODE_CONFIGURATION_SPACE_COST_HPP_
#define ROBOTOC_MULTI_MODE_CONFIGURATION_SPACE_COST_HPP_

#include <vector>
#include <unordered_map>

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
/// @class MultiModeConfigurationSpaceCost
/// @brief Multi-mode configuration space cost. 
///
class MultiModeConfigurationSpaceCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  MultiModeConfigurationSpaceCost(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  MultiModeConfigurationSpaceCost();

  ///
  /// @brief Destructor. 
  ///
  ~MultiModeConfigurationSpaceCost();

  ///
  /// @brief Default copy constructor. 
  ///
  MultiModeConfigurationSpaceCost(
      const MultiModeConfigurationSpaceCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  MultiModeConfigurationSpaceCost& operator=(
      const MultiModeConfigurationSpaceCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MultiModeConfigurationSpaceCost(
      MultiModeConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MultiModeConfigurationSpaceCost& operator=(
      MultiModeConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Sets the reference configuration q for the specified contact mode. 
  /// @param[in] q_ref Reference configuration q. Size must be Robot::dimq().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_q_ref(const Eigen::VectorXd& q_ref, const int contact_mode_id);

  ///
  /// @brief Sets the reference configuration q for the specified contact modes. 
  /// @param[in] q_ref Reference configuration q. Size must be Robot::dimq().
  /// @param[in] contact_mode_id Contact mode ids.
  ///
  void set_q_ref(const Eigen::VectorXd& q_ref, 
                 const std::vector<int>& contact_mode_id);

  ///
  /// @brief Sets the reference velocity v for the specified contact mode. 
  /// @param[in] v_ref Reference velocity v. Size must be Robot::dimv().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_v_ref(const Eigen::VectorXd& v_ref, const int contact_mode_id);

  ///
  /// @brief Sets the reference velocity v for the specified contact modes. 
  /// @param[in] v_ref Reference velocity v. Size must be Robot::dimv().
  /// @param[in] contact_mode_id Contact mode ids.
  ///
  void set_v_ref(const Eigen::VectorXd& v_ref, 
                 const std::vector<int>& contact_mode_id);

  ///
  /// @brief Sets the reference control input torques u for the specified 
  /// contact mode.
  /// @param[in] u_ref Reference control input torques u. Size must be 
  /// Robot::dimu().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_u_ref(const Eigen::VectorXd& u_ref, const int contact_mode_id);

  ///
  /// @brief Sets the reference control input torque u for the specified 
  /// contact modes. 
  /// @param[in] u_ref Reference control input torque u. Size must be 
  /// Robot::dimu().
  /// @param[in] contact_mode_id Contact mode ids.
  ///
  void set_u_ref(const Eigen::VectorXd& u_ref, 
                 const std::vector<int>& contact_mode_id);

  ///
  /// @brief Sets the weight vector on the configuration q. 
  /// @param[in] q_weight Weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_q_weight(const Eigen::VectorXd& q_weight, const int contact_mode_id);

  ///
  /// @brief Sets the weight vector on the configuration q. 
  /// @param[in] q_weight Weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_ids Contact mode ids.
  ///
  void set_q_weight(const Eigen::VectorXd& q_weight, 
                    const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight on the velocity v. 
  /// @param[in] v_weight Weight vector on the velocity v. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_v_weight(const Eigen::VectorXd& v_weight, const int contact_mode_id);

  ///
  /// @brief Sets the weight on the velocity v. 
  /// @param[in] v_weight Weight vector on the velocity v. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_ids Contact mode ids.
  ///
  void set_v_weight(const Eigen::VectorXd& v_weight, 
                    const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight on the acceleration a. 
  /// @param[in] a_weight Weight vector on the acceleration a. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_a_weight(const Eigen::VectorXd& a_weight, const int contact_mode_id);

  ///
  /// @brief Sets the weight on the acceleration a. 
  /// @param[in] a_weight Weight vector on the acceleration a. 
  /// Size must be Robot::dimv().
  /// @param[in] contact_mode_ids Contact mode ids.
  ///
  void set_a_weight(const Eigen::VectorXd& a_weight, 
                    const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight on the control input torques u. 
  /// @param[in] u_weight Weight vector on the control input torques u. 
  /// Size must be Robot::dimu().
  /// @param[in] contact_mode_id Contact mode id.
  ///
  void set_u_weight(const Eigen::VectorXd& u_weight, const int contact_mode_id);

  ///
  /// @brief Sets the weight on the control input torques u. 
  /// @param[in] u_weight Weight vector on the control input torques u. 
  /// Size must be Robot::dimu().
  /// @param[in] contact_mode_ids Contact mode ids.
  ///
  void set_u_weight(const Eigen::VectorXd& u_weight, 
                    const std::vector<int>& contact_mode_ids);

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
                          const GridInfo& grid_info, 
                          const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const GridInfo& grid_info, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const GridInfo& grid_info, 
                               const SplitSolution& s, 
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
  std::unordered_map<int, Eigen::VectorXd> q_ref_, v_ref_, u_ref_, q_weight_, 
                                           v_weight_, a_weight_, u_weight_;
  Eigen::VectorXd qf_weight_, vf_weight_, qi_weight_, vi_weight_, dvi_weight_;
};

} // namespace robotoc


#endif // ROBOTOC_MULTI_MODE_CONFIGURATION_SPACE_COST_HPP_ 