#ifndef ROBOTOC_LOCAL_CONTACT_FORCE_COST_HPP_
#define ROBOTOC_LOCAL_CONTACT_FORCE_COST_HPP_

#include <vector>

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
/// @class LocalContactForceCost
/// @brief Cost on the contact forces expressed in the local frames.
///
class LocalContactForceCost final : public CostFunctionComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  LocalContactForceCost(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  LocalContactForceCost();

  ///
  /// @brief Destructor. 
  ///
  ~LocalContactForceCost();

  ///
  /// @brief Default copy constructor. 
  ///
  LocalContactForceCost(const LocalContactForceCost&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  LocalContactForceCost& operator=(const LocalContactForceCost&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  LocalContactForceCost(LocalContactForceCost&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  LocalContactForceCost& operator=(LocalContactForceCost&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_COST_FUNCTION_COMPONENT(LocalContactForceCost)

  ///
  /// @brief Sets the reference contact forces expressed in the local frames. 
  /// @param[in] f_ref Reference contact forces expressed in the local frames. 
  /// Size must be Robot::maxNumContacts().
  ///
  void set_f_ref(const std::vector<Eigen::Vector3d>& f_ref);

  ///
  /// @brief Sets the weight vectors on the contact forces. 
  /// @param[in] f_weight Weight vectors on the contact forces. 
  /// Size must be Robot::maxNumContacts().
  ///
  void set_f_weight(const std::vector<Eigen::Vector3d>& f_weight);

  ///
  /// @brief Sets the reference impulse forces expressed in the local frames. 
  /// @param[in] fi_ref Reference impulse forces expressed in the local frames. 
  /// Size must be Robot::maxNumContacts().
  ///
  void set_fi_ref(const std::vector<Eigen::Vector3d>& fi_ref);

  ///
  /// @brief Sets the weight vectors on the impulse forces. 
  /// @param[in] fi_weight Weight vectors on the impulse forces. 
  /// Size must be Robot::maxNumContacts().
  ///
  void set_fi_weight(const std::vector<Eigen::Vector3d>& fi_weight);

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
  int max_num_contacts_, max_dimf_;
  std::vector<ContactType> contact_types_;
  std::vector<Eigen::Vector3d> f_ref_, f_weight_, fi_ref_, fi_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_LOCAL_CONTACT_FORCE_COST_HPP_ 