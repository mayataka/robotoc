#ifndef ROBOTOC_CONFIGURATION_SPACE_COST_HPP_
#define ROBOTOC_CONFIGURATION_SPACE_COST_HPP_

#include <stdexcept>

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
#include "robotoc/cost/configuration_space_ref_base.hpp"


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
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] ref Reference configuraton.
  ///
  ConfigurationSpaceCost(const Robot& robot, 
                         const std::shared_ptr<ConfigurationSpaceRefBase>& ref);

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

  std::shared_ptr<CostFunctionComponentBase> clone() const override;

  ///
  /// @brief Checks if this use the non-const reference 
  /// (ConfigurationSpaceRefBase) or not. 
  ///
  bool use_nonconst_ref() const { return use_nonconst_ref_; }

  ///
  /// @brief Sets the reference configuration. 
  /// @param[in] ref Reference configuraton.
  ///
  void set_ref(const std::shared_ptr<ConfigurationSpaceRefBase>& ref);

  ///
  /// @brief Gets the reference configuration. 
  ///
  const std::shared_ptr<ConfigurationSpaceRefBase>& get_ref() const { 
    if (!ref_) {
      throw std::runtime_error("ConfigurationSpaceRefBase is nullptr!");
    }
    return ref_; 
  }

  ///
  /// @brief Sets the const reference configuration q. 
  /// @param[in] q_ref Reference configuration q. Size must be Robot::dimq().
  ///
  void set_q_ref(const Eigen::VectorXd& q_ref);

  ///
  /// @brief Gets the const reference configuration q. 
  ///
  const Eigen::VectorXd& get_q_ref() const { return q_ref_; }

  ///
  /// @brief Sets the const reference velocity v. 
  /// @param[in] v_ref Reference velocity v. Size must be Robot::dimv().
  ///
  void set_v_ref(const Eigen::VectorXd& v_ref);

  ///
  /// @brief Gets the const reference velocity v. 
  ///
  const Eigen::VectorXd& get_v_ref() const { return v_ref_; }

  ///
  /// @brief Sets the const reference control input torques u. 
  /// @param[in] u_ref Reference control input torques u. Size must be 
  /// Robot::dimu().
  ///
  void set_u_ref(const Eigen::VectorXd& u_ref);

  ///
  /// @brief Gets the const reference control input torques u. 
  ///
  const Eigen::VectorXd& get_u_ref() const { return u_ref_; }

  ///
  /// @brief Sets the weight vector on the configuration q. 
  /// @param[in] q_weight Weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  ///
  void set_q_weight(const Eigen::VectorXd& q_weight);

  ///
  /// @brief Gets the weight vector on the configuration q. 
  ///
  const Eigen::VectorXd& get_q_weight() const { return q_weight_; }

  ///
  /// @brief Sets the weight on the velocity v. 
  /// @param[in] v_weight Weight vector on the velocity v. 
  /// Size must be Robot::dimv().
  ///
  void set_v_weight(const Eigen::VectorXd& v_weight);

  ///
  /// @brief Gets the weight on the velocity v. 
  ///
  const Eigen::VectorXd& get_v_weight() const { return v_weight_; }

  ///
  /// @brief Sets the weight on the acceleration a. 
  /// @param[in] a_weight Weight vector on the acceleration a. 
  /// Size must be Robot::dimv().
  ///
  void set_a_weight(const Eigen::VectorXd& a_weight);

  ///
  /// @brief Gets the weight on the acceleration a. 
  ///
  const Eigen::VectorXd& get_a_weight() const { return a_weight_; }

  ///
  /// @brief Sets the weight on the control input torques u. 
  /// @param[in] u_weight Weight vector on the control input torques u. 
  /// Size must be Robot::dimu().
  ///
  void set_u_weight(const Eigen::VectorXd& u_weight);

  ///
  /// @brief Gets the weight on the control input torques u. 
  ///
  const Eigen::VectorXd& get_u_weight() const { return u_weight_; }

  ///
  /// @brief Sets the weight vector on the configuration q at the terminal stage. 
  /// @param[in] q_weight_terminal Weight vector on the configuration q at the terminal 
  /// stage. Size must be Robot::dimv().
  ///
  void set_q_weight_terminal(const Eigen::VectorXd& q_weight_terminal);

  ///
  /// @brief Gets the weight vector on the configuration q at the terminal stage. 
  ///
  const Eigen::VectorXd& get_q_weight_terminal() const { return q_weight_terminal_; }

  ///
  /// @brief Sets the weight vector on the velocity v at the terminal stage. 
  /// @param[in] v_weight_terminal Weight vector on the velocity v at the terminal 
  /// stage. Size must be Robot::dimv().
  ///
  void set_v_weight_terminal(const Eigen::VectorXd& v_weight_terminal);

  ///
  /// @brief Gets the weight vector on the velocity v at the terminal stage. 
  ///
  const Eigen::VectorXd& get_v_weight_terminal() const { return v_weight_terminal_; }

  ///
  /// @brief Sets the weight vector on the configuration q at impulse stages. 
  /// @param[in] q_weight_impulse Weight vector on the configuration q at impulse  
  /// stages. Size must be Robot::dimv().
  ///
  void set_q_weight_impulse(const Eigen::VectorXd& q_weight_impulse);

  ///
  /// @brief Gets the weight vector on the configuration q at impulse stages. 
  ///
  const Eigen::VectorXd& get_q_weight_impulse() const { return q_weight_impulse_; }

  ///
  /// @brief Sets the weight vector on the velocity v at the impulse stages. 
  /// @param[in] v_weight_impulse Weight vector on the velocity v at the impulse  
  /// stages. Size must be Robot::dimv().
  ///
  void set_v_weight_impulse(const Eigen::VectorXd& v_weight_impulse);

  ///
  /// @brief Gets the weight vector on the velocity v at the impulse stages. 
  ///
  const Eigen::VectorXd& get_v_weight_impulse() const { return v_weight_impulse_; }

  ///
  /// @brief Sets the weight vector on the impulse change in the velocity dv at 
  /// the impulse stages. 
  /// @param[in] dv_weight_impulse Weight vector on the impulse change in the velocity
  /// the impulse stages. Size must be Robot::dimv().
  ///
  void set_dv_weight_impulse(const Eigen::VectorXd& dv_weight_impulse);

  ///
  /// @brief Gets the weight vector on the impulse change in the velocity dv at 
  /// the impulse stages. 
  ///
  const Eigen::VectorXd& get_dv_weight_impulse() const { return dv_weight_impulse_; }

  ///
  /// @brief Sets (clones) reference and weight parameters from other. 
  /// @param[in] other other config cost. 
  ///
  void set_from_other(const std::shared_ptr<ConfigurationSpaceCost>& other);

  ///
  /// @brief Evaluate if the cost on the configuration q is active for given 
  /// grid_info. 
  /// @param[in] grid_info Grid info.
  /// @return Cost status (if the cost is active or not).
  ///
  bool isCostConfigActive(const GridInfo& grid_info) const {
    if (use_nonconst_ref_) {
      return ref_->isActive(grid_info);
    }
    else {
      return true;
    }
  }

  ///
  /// @brief Evaluate the difference between the configuration and the reference
  /// configuration. 
  /// @param[in] robot Robot model.
  /// @param[in, out] data Cost funciton data.
  /// @param[in] grid_info Grid info
  /// @param[in] q Current configuration 
  ///
  void evalConfigDiff(const Robot& robot, CostFunctionData& data, 
                      const GridInfo& grid_info, 
                      const Eigen::VectorXd& q) const {
    if (use_nonconst_ref_) {
      if (ref_->isActive(grid_info)) {
        ref_->updateRef(robot, grid_info, data.q_ref);
        robot.subtractConfiguration(q, data.q_ref, data.qdiff);
      }
    }
    else {
      robot.subtractConfiguration(q, q_ref_, data.qdiff);
    }
  }

  ///
  /// @brief Evaluate the Jacobian of the difference between the configuration 
  /// and the reference configuration. 
  /// @param[in] robot Robot model.
  /// @param[in, out] data Cost funciton data.
  /// @param[in] grid_info Grid info
  /// @param[in] q Current configuration 
  ///
  void evalConfigDiffJac(const Robot& robot, CostFunctionData& data, 
                         const GridInfo& grid_info, 
                         const Eigen::VectorXd& q) const {
    if (use_nonconst_ref_) {
      if (ref_->isActive(grid_info)) {
        robot.dSubtractConfiguration_dqf(q, data.q_ref, data.J_qdiff);
      }
    }
    else {
      robot.dSubtractConfiguration_dqf(q, q_ref_, data.J_qdiff);
    }
  }

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
  Eigen::VectorXd q_ref_, v_ref_, u_ref_, 
                  q_weight_, v_weight_, a_weight_, u_weight_,
                  q_weight_terminal_, v_weight_terminal_, 
                  q_weight_impulse_, v_weight_impulse_, dv_weight_impulse_;
  std::shared_ptr<ConfigurationSpaceRefBase> ref_;
  bool use_nonconst_ref_, 
       enable_q_cost_, enable_v_cost_, enable_a_cost_, enable_u_cost_,
       enable_q_cost_terminal_, enable_v_cost_terminal_,
       enable_q_cost_impulse_, enable_v_cost_impulse_, enable_dv_cost_impulse_;
};

} // namespace robotoc


#endif // ROBOTOC_CONFIGURATION_SPACE_COST_HPP_ 