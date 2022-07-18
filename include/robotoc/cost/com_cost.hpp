#ifndef ROBOTOC_COM_COST_HPP_
#define ROBOTOC_COM_COST_HPP_

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
#include "robotoc/cost/com_ref_base.hpp"


namespace robotoc {

///
/// @class CoMCost
/// @brief Cost on the position of the center of mass. 
///
class CoMCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  CoMCost(const Robot& robot);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] ref Reference CoM position.
  ///
  CoMCost(const Robot& robot, const std::shared_ptr<CoMRefBase>& ref);

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] const_ref Const reference CoM position.
  ///
  CoMCost(const Robot& robot, const Eigen::Vector3d& const_ref);

  ///
  /// @brief Default constructor. 
  ///
  CoMCost();

  ///
  /// @brief Destructor. 
  ///
  ~CoMCost();

  ///
  /// @brief Default copy constructor. 
  ///
  CoMCost(const CoMCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CoMCost& operator=(const CoMCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CoMCost(CoMCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CoMCost& operator=(CoMCost&&) noexcept = default;

  std::shared_ptr<CostFunctionComponentBase> clone() const override;

  ///
  /// @brief Checks if this use the non-const reference 
  /// (CoMRefBase) or not. 
  ///
  bool use_nonconst_ref() const { return use_nonconst_ref_; }

  ///
  /// @brief Sets the reference CoM position. 
  /// @param[in] ref Reference CoM position.
  ///
  void set_ref(const std::shared_ptr<CoMRefBase>& ref);

  ///
  /// @brief Gets the reference CoM position.
  ///
  const std::shared_ptr<CoMRefBase>& get_ref() const { return ref_; }

  ///
  /// @brief Sets the const reference CoM position. 
  /// @param[in] const_ref Const reference CoM position.
  ///
  void set_const_ref(const Eigen::Vector3d& const_ref);

  ///
  /// @brief Gets the const reference CoM position. 
  ///
  const Eigen::Vector3d& get_const_ref() const { return const_ref_; } 

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] weight Weight vector on the CoM position error. 
  ///
  void set_weight(const Eigen::Vector3d& weight);

  ///
  /// @brief Gets the weight vector. 
  ///
  const Eigen::Vector3d& get_weight() const { return weight_; } 

  ///
  /// @brief Sets the weight vector for the terminal stage. 
  /// @param[in] weight_terminal Weight vector on the CoM position error 
  /// at the terminal stage. 
  ///
  void set_weight_terminal(const Eigen::Vector3d& weight_terminal);

  ///
  /// @brief Gets the weight vector for the terminal stage. 
  ///
  const Eigen::Vector3d& get_weight_terminal() const { return weight_terminal_; } 

  ///
  /// @brief Sets the weight vector for the impulse stage. 
  /// @param[in] weight_impulse Weight vector on the CoM position error 
  /// at the impulse stage. 
  ///
  void set_weight_impulse(const Eigen::Vector3d& weight_impulse);

  ///
  /// @brief Gets the weight vector for the impulse stage. 
  ///
  const Eigen::Vector3d& get_weight_impulse() const { return weight_impulse_; } 

  ///
  /// @brief Sets (clones) reference and weight parameters from other. 
  /// @param[in] other other CoM cost. 
  ///
  void set_from_other(const std::shared_ptr<CoMCost>& other);

  ///
  /// @brief Evaluate if the cost is active for given grid_info. 
  /// @param[in] grid_info Grid info.
  /// @return Cost status (if the cost is active or not).
  ///
  bool isCostActive(const GridInfo& grid_info) const {
    if (use_nonconst_ref_) {
      return ref_->isActive(grid_info);
    }
    else {
      return true;
    }
  }

  ///
  /// @brief Evaluate the difference between the robot's task-space position 
  /// status and reference. 
  /// @param[in] robot Robot model.
  /// @param[in, out] data Cost funciton data.
  /// @param[in] grid_info Grid info
  ///
  void evalDiff(const Robot& robot, CostFunctionData& data, 
                const GridInfo& grid_info) const {
    if (use_nonconst_ref_) {
      if (ref_->isActive(grid_info)) {
        ref_->updateRef(grid_info, data.x3d_ref);
        data.diff_3d = robot.CoM() - data.x3d_ref;
      }
    }
    else {
      data.diff_3d = robot.CoM() - const_ref_;
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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::Vector3d const_ref_, weight_, weight_terminal_, weight_impulse_;
  std::shared_ptr<CoMRefBase> ref_;
  bool use_nonconst_ref_, enable_cost_, enable_cost_terminal_, enable_cost_impulse_;
};

} // namespace robotoc

#endif // ROBOTOC_COM_COST_HPP_ 