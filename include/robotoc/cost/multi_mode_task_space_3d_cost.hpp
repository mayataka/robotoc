#ifndef ROBOTOC_MULTI_MODE_TASK_SPACE_3D_COST_HPP_
#define ROBOTOC_MULTI_MODE_TASK_SPACE_3D_COST_HPP_

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
/// @class MultiModeTaskSpace3DCost
/// @brief Multi-mode cost on the task space position. 
///
class MultiModeTaskSpace3DCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  ///
  MultiModeTaskSpace3DCost(const Robot& robot, const int frame_id);

  ///
  /// @brief Default constructor. 
  ///
  MultiModeTaskSpace3DCost();

  ///
  /// @brief Destructor. 
  ///
  ~MultiModeTaskSpace3DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  MultiModeTaskSpace3DCost(const MultiModeTaskSpace3DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  MultiModeTaskSpace3DCost& operator=(const MultiModeTaskSpace3DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MultiModeTaskSpace3DCost(MultiModeTaskSpace3DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MultiModeTaskSpace3DCost& operator=(MultiModeTaskSpace3DCost&&) noexcept = default;

  ///
  /// @brief Sets the reference position. 
  /// @param[in] x3d_ref Reference position.
  /// @param[in] contact_mode_id Contact mode id. Default is 0.
  ///
  void set_x3d_ref(const Eigen::Vector3d& x3d_ref, const int contact_mode_id=0);

  ///
  /// @brief Sets the reference position. 
  /// @param[in] x3d_ref Reference position.
  /// @param[in] contact_mode_ids Contact mode ids. 
  ///
  void set_x3d_ref(const Eigen::Vector3d& x3d_ref, 
                   const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] x3d_weight Weight vector on the position error. 
  /// @param[in] contact_mode_id Contact mode id. Default is 0.
  ///
  void set_x3d_weight(const Eigen::Vector3d& x3d_weight, 
                      const int contact_mode_id=0);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] x3d_weight Weight vector on the position error. 
  /// @param[in] contact_mode_ids Contact mode ids. 
  ///
  void set_x3d_weight(const Eigen::Vector3d& x3d_weight, 
                      const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight vector at the terminal stage. 
  /// @param[in] x3df_weight Weight vector on the position error at the
  /// terminal stage. 
  ///
  void set_x3df_weight(const Eigen::Vector3d& x3df_weight);

  ///
  /// @brief Sets the weight vector at the impulse stage. 
  /// @param[in] x3di_weight Weight vector on the position error at the 
  /// impulse stage. 
  ///
  void set_x3di_weight(const Eigen::Vector3d& x3di_weight);

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
  std::unordered_map<int, Eigen::Vector3d> x3d_ref_, x3d_weight_;
  Eigen::Vector3d x3df_weight_, x3di_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_MULTI_MODE_TASK_SPACE_3D_COST_HPP_ 