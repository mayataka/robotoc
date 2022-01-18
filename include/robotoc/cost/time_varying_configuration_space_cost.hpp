#ifndef ROBOTOC_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_
#define ROBOTOC_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_

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
/// @class TimeVaryingConfigurationRefBase 
/// @brief Base class of time-varying reference of the configuration. 
///
class TimeVaryingConfigurationRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingConfigurationRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TimeVaryingConfigurationRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingConfigurationRefBase(
      const TimeVaryingConfigurationRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingConfigurationRefBase& operator=(
      const TimeVaryingConfigurationRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingConfigurationRefBase(
      TimeVaryingConfigurationRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingConfigurationRefBase& operator=(
      TimeVaryingConfigurationRefBase&&) noexcept = default;

  ///
  /// @brief Computes the time-varying reference configuration. 
  /// @param[in] robot Robot model.
  /// @param[in] t Time.
  /// @param[in] q_ref Reference position. Size is Robot::dimv().
  ///
  virtual void update_q_ref(const Robot& robot, const double t, 
                            Eigen::VectorXd& q_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] t Time.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const double t) const = 0;

};


///
/// @class TimeVaryingConfigurationSpaceCost
/// @brief Cost on the time-varying configuration. 
///
class TimeVaryingConfigurationSpaceCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] ref Shared ptr to the reference configuration.
  ///
  TimeVaryingConfigurationSpaceCost(
      const Robot& robot, 
      const std::shared_ptr<TimeVaryingConfigurationRefBase>& ref);

  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingConfigurationSpaceCost();

  ///
  /// @brief Destructor. 
  ///
  ~TimeVaryingConfigurationSpaceCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingConfigurationSpaceCost(
      const TimeVaryingConfigurationSpaceCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingConfigurationSpaceCost& operator=(
      const TimeVaryingConfigurationSpaceCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingConfigurationSpaceCost(
      TimeVaryingConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingConfigurationSpaceCost& operator=(
      TimeVaryingConfigurationSpaceCost&&) noexcept = default;

  ///
  /// @brief Sets the time-varying reference configuration. 
  /// @param[in] q_ref Shared ptr to the time-varying reference position.
  ///
  void set_q_ref(const std::shared_ptr<TimeVaryingConfigurationRefBase>& q_ref);

  ///
  /// @brief Sets the weight vector on the configuration q. 
  /// @param[in] q_weight Weight vector on the configuration q. 
  /// Size must be Robot::dimv().
  ///
  void set_q_weight(const Eigen::VectorXd& q_weight);

  ///
  /// @brief Sets the weight vector on the configuration q at the terminal stage. 
  /// @param[in] qf_weight Weight vector on the configuration q at the terminal 
  /// stage. Size must be Robot::dimv().
  ///
  void set_qf_weight(const Eigen::VectorXd& qf_weight);

  ///
  /// @brief Sets the weight vector on the configuration q at the impulse stages. 
  /// @param[in] qi_weight Weight vector on the configuration q at the impulse 
  /// stages. Size must be Robot::dimv().
  ///
  void set_qi_weight(const Eigen::VectorXd& qi_weight);

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
                                  ImpulseSplitKKTResidual& kkt_residual) const;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const GridInfo& grid_info,
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

private:
  int dimq_, dimv_;
  std::shared_ptr<TimeVaryingConfigurationRefBase> q_ref_;
  Eigen::VectorXd q_weight_, qf_weight_, qi_weight_;

};

} // namespace robotoc

#endif // ROBOTOC_TIME_VARYING_CONFIGURATION_SPACE_COST_HPP_ 