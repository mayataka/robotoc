#ifndef IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_
#define IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_base.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class TimeVaryingTaskSpace3DRefBase
/// @brief Base class of time-varying reference of task space position. 
///
class TimeVaryingTaskSpace3DRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingTaskSpace3DRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TimeVaryingTaskSpace3DRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingTaskSpace3DRefBase(const TimeVaryingTaskSpace3DRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingTaskSpace3DRefBase& operator=(
      const TimeVaryingTaskSpace3DRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingTaskSpace3DRefBase(
      TimeVaryingTaskSpace3DRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingTaskSpace3DRefBase& operator=(
      TimeVaryingTaskSpace3DRefBase&&) noexcept = default;

  ///
  /// @brief Computes the time-varying reference position. 
  /// @param[in] t Time.
  /// @param[in] q_3d_ref Reference position. Size is 3.
  ///
  virtual void update_q_3d_ref(const double t, 
                               Eigen::VectorXd& q_3d_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] t Time.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const double t) const = 0;

};


///
/// @class TimeVaryingTaskSpace3DCost 
/// @brief Cost on the time-varying task space position. 
///
class TimeVaryingTaskSpace3DCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  /// @param[in] ref Shared ptr to the reference position.
  ///
  TimeVaryingTaskSpace3DCost(
      const Robot& robot, const int frame_id, 
      const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& ref);

  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingTaskSpace3DCost();

  ///
  /// @brief Destructor. 
  ///
  ~TimeVaryingTaskSpace3DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingTaskSpace3DCost(const TimeVaryingTaskSpace3DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingTaskSpace3DCost& operator=(
      const TimeVaryingTaskSpace3DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingTaskSpace3DCost(TimeVaryingTaskSpace3DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingTaskSpace3DCost& operator=(
      TimeVaryingTaskSpace3DCost&&) noexcept = default;

  ///
  /// @brief Sets the time-varying reference position. 
  /// @param[in] ref Shared ptr to the time-varying reference position.
  ///
  void set_ref(const std::shared_ptr<TimeVaryingTaskSpace3DRefBase>& ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] q_3d_weight Weight vector on the position error. 
  ///
  void set_q_weight(const Eigen::Vector3d& q_3d_weight);

  ///
  /// @brief Sets the terminal weight vector. 
  /// @param[in] qf_3d_weight Terminal weight vector on the position error. 
  ///
  void set_qf_weight(const Eigen::Vector3d& qf_3d_weight);

  ///
  /// @brief Sets the weight vector at impulse. 
  /// @param[in] qi_3d_weight Weight vector on the position error at impulse. 
  ///
  void set_qi_weight(const Eigen::Vector3d& qi_3d_weight);

  bool useKinematics() const override;

  double computeStageCost(Robot& robot, CostFunctionData& data, const double t, 
                          const double dt, 
                          const SplitSolution& s) const override;

  double computeTerminalCost(Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const SplitSolution& s) const override;

  double computeImpulseCost(Robot& robot, CostFunctionData& data, 
                            const double t, 
                            const ImpulseSplitSolution& s) const override;

  void computeStageCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, const double dt, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeTerminalCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const SplitSolution& s, SplitKKTResidual& kkt_residual) const override;

  void computeImpulseCostDerivatives(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const;

  void computeStageCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const double dt, 
                               const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  void computeTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                                  const double t, const SplitSolution& s, 
                                  SplitKKTMatrix& kkt_matrix) const override;

  void computeImpulseCostHessian(
      Robot& robot, CostFunctionData& data, const double t, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  std::shared_ptr<TimeVaryingTaskSpace3DRefBase> ref_;
  Eigen::Vector3d q_3d_weight_, qf_3d_weight_, qi_3d_weight_;

};

} // namespace idocp


#endif // IDOCP_TIME_VARYING_TASK_SPACE_3D_COST_HPP_ 