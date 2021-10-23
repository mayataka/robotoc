#ifndef ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_
#define ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


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
  /// @param[in] t Time.
  /// @param[in] SE3_ref Reference placement.
  ///
  virtual void update_SE3_ref(const double t, SE3& SE3_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] t Time.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const double t) const = 0;

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
  /// @param[in] ref Shared ptr time-varying reference placement.
  ///
  void set_ref(const std::shared_ptr<TimeVaryingTaskSpace6DRefBase>& ref);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] position_weight Weight vector on the position error. 
  /// @param[in] rotation_weight Weight vector on the rotation error. 
  ///
  void set_q_weight(const Eigen::Vector3d& position_weight, 
                       const Eigen::Vector3d& rotation_weight);

  ///
  /// @brief Sets the terminal weight vectors. 
  /// @param[in] position_weight Temrinal weight vector on the position error. 
  /// @param[in] rotation_weight Temrinal weight vector on the rotation error. 
  ///
  void set_qf_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  ///
  /// @brief Sets the weight vectors at impulse. 
  /// @param[in] position_weight Weight vector on the position error at impulse. 
  /// @param[in] rotation_weight Weight vector on the rotation error at impulse. 
  ///
  void set_qi_weight(const Eigen::Vector3d& position_weight, 
                        const Eigen::Vector3d& rotation_weight);

  bool useKinematics() const override;

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override;

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override;

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const override;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override;

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const;

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  std::shared_ptr<TimeVaryingTaskSpace6DRefBase> ref_;
  Eigen::VectorXd q_6d_weight_, qf_6d_weight_, qi_6d_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_TIME_VARYING_TASK_SPACE_6D_COST_HPP_ 