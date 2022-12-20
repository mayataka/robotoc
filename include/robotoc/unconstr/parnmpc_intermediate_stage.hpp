#ifndef ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_
#define ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/unconstr/unconstr_ocp_data.hpp"
#include "robotoc/ocp/grid_info.hpp"


namespace robotoc {

///
/// @class ParNMPCIntermediateStage
/// @brief The intermediate stage of ParNMPC for unconstrained rigid-body systems.
///
class ParNMPCIntermediateStage {
public:
  ///
  /// @brief Constructs the intermediate stage.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  ParNMPCIntermediateStage(const Robot& robot, 
                           const std::shared_ptr<CostFunction>& cost,
                           const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  ParNMPCIntermediateStage();

  ///
  /// @brief Default destructor. 
  ///
  ~ParNMPCIntermediateStage() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCIntermediateStage(const ParNMPCIntermediateStage&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCIntermediateStage& operator=(const ParNMPCIntermediateStage&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCIntermediateStage(ParNMPCIntermediateStage&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCIntermediateStage& operator=(ParNMPCIntermediateStage&&) noexcept = default;

  ///
  /// @brief Creates the data.
  /// @param[in] robot Robot model. 
  ///
  UnconstrOCPData createData(const Robot& robot) const;

  ///
  /// @brief Checks whether the solution is feasible w.r.t. the inequality constraints.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  bool isFeasible(Robot& robot, const GridInfo& grid_info,
                  const SplitSolution& s, UnconstrOCPData& data) const;

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  void initConstraints(Robot& robot, const GridInfo& grid_info, 
                       const SplitSolution& s, UnconstrOCPData& data) const;

  ///
  /// @brief Computes the stage cost and constraint violation of this stage.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info. 
  /// @param[in] q_prev Configuration at the previous stage. 
  /// @param[in] v_prev Configuration at the previous stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this stage.
  ///
  void evalOCP(Robot& robot, const GridInfo& grid_info, 
               const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
               const SplitSolution& s, UnconstrOCPData& data, 
               SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the KKT matrix and residual of this stage.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info. 
  /// @param[in] q_prev Configuration at the previous stage. 
  /// @param[in] v_prev Configuration at the previous stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this stage.
  /// @param[in, out] kkt_residual Split KKT residual of this stage.
  ///
  void evalKKT(Robot& robot, const GridInfo& grid_info, 
               const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
               const SplitSolution& s, const SplitSolution& s_next, 
               UnconstrOCPData& data, SplitKKTMatrix& kkt_matrix, 
               SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Expands the primal and dual variables, i.e., computes the Newton 
  /// direction of the condensed variables of this stage.
  /// @param[in] grid_info Grid info. 
  /// @param[in] kkt_matrix Split KKT matrix of this stage.
  /// @param[in] kkt_residual Split KKT residual of this stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] d Split direction of this stage.
  ///
  void expandPrimalAndDual(const GridInfo& grid_info,
                           const SplitKKTMatrix& kkt_matrix, 
                           const SplitKKTResidual& kkt_residual, 
                           UnconstrOCPData& data, SplitDirection& d) const;

  ///
  /// @brief Computes the maximum primal step size.
  /// @param[in] data Data of this stage. 
  /// @return Maximum primal step size.
  ///
  double maxPrimalStepSize(const UnconstrOCPData& data) const;

  ///
  /// @brief Computes the maximum dual size.
  /// @param[in] data Data of this stage. 
  /// @return Maximum dual step size.
  ///
  double maxDualStepSize(const UnconstrOCPData& data) const;

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s, 
                    UnconstrOCPData& data) const;

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size. 
  /// @param[in, out] data Data of this stage. 
  ///
  void updateDual(const double dual_step_size, UnconstrOCPData& data) const;

private:
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  ContactStatus contact_status_;
};

} // namespace robotoc

#endif // ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_ 