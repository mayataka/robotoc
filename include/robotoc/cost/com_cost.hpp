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

  ///
  /// @brief Sets the reference position of the center of mass. 
  /// @param[in] com_ref Reference position of the center of mass.
  ///
  void set_com_ref(const Eigen::Vector3d& com_ref);

  ///
  /// @brief Sets the weight vector. 
  /// @param[in] com_weight Weight vector on the com position error. 
  ///
  void set_com_weight(const Eigen::Vector3d& com_weight);

  ///
  /// @brief Sets the weight vector at the terminal stage. 
  /// @param[in] comf_weight Weight vector on the com position error at the 
  /// terminal stage. 
  ///
  void set_comf_weight(const Eigen::Vector3d& comf_weight);

  ///
  /// @brief Sets the weight vector at the impulse stage. 
  /// @param[in] comi_weight Weight vector on the com position error at the 
  /// impulse stage. 
  ///
  void set_comi_weight(const Eigen::Vector3d& comi_weight);

  bool useKinematics() const override;

  double evalStageCost(Robot& robot, const ContactStatus& contact_status, 
                       CostFunctionData& data, const int time_stage_in_phase, 
                       const double t, const double dt, 
                       const SplitSolution& s) const override;

  void evalStageCostDerivatives(Robot& robot, const ContactStatus& contact_status, 
                                CostFunctionData& data, const int time_stage_in_phase, 
                                const double t, const double dt, 
                                const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual) const override;

  void evalStageCostHessian(Robot& robot, const ContactStatus& contact_status, 
                            CostFunctionData& data, const int time_stage_in_phase, 
                            const double t, const double dt, 
                            const SplitSolution& s, 
                            SplitKKTMatrix& kkt_matrix) const override;

  double evalTerminalCost(Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s) const override;

  void evalTerminalCostDerivatives(Robot& robot, CostFunctionData& data, 
                                   const double t, const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const override;

  void evalTerminalCostHessian(Robot& robot, CostFunctionData& data, 
                               const double t, const SplitSolution& s, 
                               SplitKKTMatrix& kkt_matrix) const override;

  double evalImpulseCost(Robot& robot, const ImpulseStatus& impulse_status, 
                         CostFunctionData& data, const double t, 
                         const ImpulseSplitSolution& s) const override;

  void evalImpulseCostDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                                  CostFunctionData& data, const double t, 
                                  const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTResidual& kkt_residual) const;

  void evalImpulseCostHessian(Robot& robot, const ImpulseStatus& impulse_status, 
                              CostFunctionData& data, const double t, 
                              const ImpulseSplitSolution& s, 
                              ImpulseSplitKKTMatrix& kkt_matrix) const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::Vector3d com_ref_, com_weight_, comf_weight_, comi_weight_;

};

} // namespace robotoc


#endif // ROBOTOC_COM_COST_HPP_ 