#ifndef ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_
#define ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

///
/// @class TimeVaryingTaskSpace6DRefBase
/// @brief Base class of time-varying reference of task space placement. 
///
class TimeVaryingTaskSpace6DRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingTaskSpace6DRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TimeVaryingTaskSpace6DRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingTaskSpace6DRefBase(const TimeVaryingTaskSpace6DRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingTaskSpace6DRefBase& operator=(
      const TimeVaryingTaskSpace6DRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingTaskSpace6DRefBase(
      TimeVaryingTaskSpace6DRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingTaskSpace6DRefBase& operator=(
      TimeVaryingTaskSpace6DRefBase&&) noexcept = default;

  ///
  /// @brief Computes the time-varying reference placement. 
  /// @param[in] grid_info Grid info.
  /// @param[in] x6d_ref Reference placement.
  ///
  virtual void update_x6d_ref(const GridInfo& grid_info, 
                              SE3& x6d_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] grid_info Grid info.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const GridInfo& grid_info) const = 0;
};


///
/// @class TimeVaryingTaskSpace6DCost 
/// @brief Cost on the time-varying task space placement. 
///
class TimeVaryingTaskSpace6DCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] ref Shared ptr to the reference placement.
  ///
  TimeVaryingTaskSpace6DCost(
      const Robot& robot, const int frame_id, 
      const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref);

  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingTaskSpace6DCost();

  ///
  /// @brief Destructor. 
  ///
  ~TimeVaryingTaskSpace6DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingTaskSpace6DCost(const TimeVaryingTaskSpace6DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingTaskSpace6DCost& operator=(
      const TimeVaryingTaskSpace6DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingTaskSpace6DCost(TimeVaryingTaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingTaskSpace6DCost& operator=(
      TimeVaryingTaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Sets the time-varying reference placement. 
  /// @param[in] x6d_ref Shared ptr time-varying reference placement.
  ///
  void set_x6d_ref(
      const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& x6d_ref);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] trans_weight Weight vector on the position error. 
  /// @param[in] rot_weight Weight vector on the rotation error. 
  ///
  void set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                      const Eigen::Vector3d& rot_weight);

  ///
  /// @brief Sets the weight vectors at the terminal stage. 
  /// @param[in] trans_weight Weight vector on the position error at the 
  /// terminal stage. 
  /// @param[in] rot_weight Weight vector on the rotation error at the terminal 
  /// stage. 
  ///
  void set_x6df_weight(const Eigen::Vector3d& trans_weight, 
                       const Eigen::Vector3d& rot_weight);

  ///
  /// @brief Sets the weight vectors at the impulse stages. 
  /// @param[in] trans_weight Weight vector on the position error at the 
  /// impulse stages. 
  /// @param[in] rot_weight Weight vector on the rotation error at the 
  /// impulse stages. 
  ///
  void set_x6di_weight(const Eigen::Vector3d& trans_weight, 
                       const Eigen::Vector3d& rot_weight);

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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  std::shared_ptr<TimeVaryingTaskSpace6DRefBase> x6d_ref_;
  Eigen::VectorXd x6d_weight_, x6df_weight_, x6di_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_ 