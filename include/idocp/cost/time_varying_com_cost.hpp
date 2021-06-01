#ifndef IDOCP_TIME_VARYING_COM_COST_HPP_
#define IDOCP_TIME_VARYING_COM_COST_HPP_

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
/// @class TimeVaryingCoMRefBase
/// @brief Base class of time-varying reference of position of the center of mass. 
///
class TimeVaryingCoMRefBase {
public:
  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingCoMRefBase() {}

  ///
  /// @brief Destructor. 
  ///
  virtual ~TimeVaryingCoMRefBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingCoMRefBase(const TimeVaryingCoMRefBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingCoMRefBase& operator=(
      const TimeVaryingCoMRefBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingCoMRefBase(
      TimeVaryingCoMRefBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingCoMRefBase& operator=(
      TimeVaryingCoMRefBase&&) noexcept = default;

  ///
  /// @brief Computes the time-varying reference position. 
  /// @param[in] t Time.
  /// @param[in] CoM_ref Reference position of the center of mass. Size is 3.
  ///
  virtual void update_CoM_ref(const double t, 
                              Eigen::VectorXd& CoM_ref) const = 0;

  ///
  /// @brief Checks wheather the cost is active or not at the specified time. 
  /// @param[in] t Time.
  /// @return true if the cost is active at time t. false if not.
  ///
  virtual bool isActive(const double t) const = 0;

};


///
/// @class TimeVaryingCoMCost 
/// @brief Cost on the time-varying position of the center of mass. 
///
class TimeVaryingCoMCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] ref Shared ptr to the reference position.
  ///
  TimeVaryingCoMCost(const Robot& robot, 
                     const std::shared_ptr<TimeVaryingCoMRefBase>& ref);

  ///
  /// @brief Default constructor. 
  ///
  TimeVaryingCoMCost();

  ///
  /// @brief Destructor. 
  ///
  ~TimeVaryingCoMCost();

  ///
  /// @brief Default copy constructor. 
  ///
  TimeVaryingCoMCost(const TimeVaryingCoMCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  TimeVaryingCoMCost& operator=(const TimeVaryingCoMCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TimeVaryingCoMCost(TimeVaryingCoMCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TimeVaryingCoMCost& operator=(TimeVaryingCoMCost&&) noexcept = default;

  ///
  /// @brief Sets the time-varying reference position. 
  /// @param[in] ref Shared ptr to the time-varying reference position.
  ///
  void set_ref(const std::shared_ptr<TimeVaryingCoMRefBase>& ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] q_weight Weight vector on the CoM position error. 
  ///
  void set_q_weight(const Eigen::Vector3d& q_weight);

  ///
  /// @brief Sets the terminal weight vector. 
  /// @param[in] qf_weight Terminal weight vector on the CoM position error. 
  ///
  void set_qf_weight(const Eigen::Vector3d& qf_weight);

  ///
  /// @brief Sets the weight vector at impulse. 
  /// @param[in] qi_weight Weight vector on the CoM position error at impulse. 
  ///
  void set_qi_weight(const Eigen::Vector3d& qi_weight);

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
  std::shared_ptr<TimeVaryingCoMRefBase> ref_;
  Eigen::Vector3d q_weight_, qf_weight_, qi_weight_;

};

} // namespace idocp


#endif // IDOCP_TIME_VARYING_COM_COST_HPP_ 