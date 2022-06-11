#ifndef ROBOTOC_MULTI_MODE_TASK_SPACE_6D_COST_HPP_
#define ROBOTOC_MULTI_MODE_TASK_SPACE_6D_COST_HPP_

#include <vector>
#include <unordered_map>

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


namespace robotoc {

///
/// @class MultiModeTaskSpace6DCost
/// @brief Cost on the task space placement (SE(3)). 
///
class MultiModeTaskSpace6DCost final : public CostFunctionComponentBase {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] frame_id Frame of interest.
  ///
  MultiModeTaskSpace6DCost(const Robot& robot, const int frame_id);

  ///
  /// @brief Default constructor. 
  ///
  MultiModeTaskSpace6DCost();

  ///
  /// @brief Destructor. 
  ///
  ~MultiModeTaskSpace6DCost();

  ///
  /// @brief Default copy constructor. 
  ///
  MultiModeTaskSpace6DCost(const MultiModeTaskSpace6DCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  MultiModeTaskSpace6DCost& operator=(const MultiModeTaskSpace6DCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MultiModeTaskSpace6DCost(MultiModeTaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MultiModeTaskSpace6DCost& operator=(MultiModeTaskSpace6DCost&&) noexcept = default;

  ///
  /// @brief Sets the reference pose. 
  /// @param[in] trans_ref Reference translation.
  /// @param[in] rot_ref Reference rotation matrix.
  /// @param[in] contact_mode_id Contact mode id. Default is 0.
  ///
  void set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                   const Eigen::Matrix3d& rot_ref,
                   const int contact_mode_id=0);

  ///
  /// @brief Sets the reference pose. 
  /// @param[in] trans_ref Reference translation.
  /// @param[in] rot_ref Reference rotation matrix.
  /// @param[in] contact_mode_ids Contact mode ids. 
  ///
  void set_x6d_ref(const Eigen::Vector3d& trans_ref, 
                   const Eigen::Matrix3d& rot_ref,
                   const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] trans_weight Weight vector on the position error. 
  /// @param[in] rot_weight Weight vector on the rotation error. 
  /// @param[in] contact_mode_id Contact mode id. Default is 0.
  ///
  void set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                      const Eigen::Vector3d& rot_weight,
                      const int contact_mode_id=0);

  ///
  /// @brief Sets the weight vectors. 
  /// @param[in] trans_weight Weight vector on the position error. 
  /// @param[in] rot_weight Weight vector on the rotation error. 
  /// @param[in] contact_mode_ids Contact mode id.
  ///
  void set_x6d_weight(const Eigen::Vector3d& trans_weight, 
                      const Eigen::Vector3d& rot_weight,
                      const std::vector<int>& contact_mode_ids);

  ///
  /// @brief Sets the weight vectors at the terminal stage. 
  /// @param[in] trans_weight Weight vector on the position error at the 
  /// terminal stage. 
  /// @param[in] rot_weight Weight vector on the rotation error at the 
  /// terminal stage.
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
                                  ImpulseSplitKKTResidual& kkt_residual) const override;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const GridInfo& grid_info,
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int frame_id_;
  std::unordered_map<int, SE3> x6d_ref_, x6d_ref_inv_;
  std::unordered_map<int, Vector6d> x6d_weight_;
  Vector6d x6df_weight_, x6di_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_MULTI_MODE_TASK_SPACE_6D_COST_HPP_ 