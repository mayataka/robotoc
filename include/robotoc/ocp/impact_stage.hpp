#ifndef ROBOTOC_IMPACT_STAGE_HPP_
#define ROBOTOC_IMPACT_STAGE_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/ocp/ocp_data.hpp"


namespace robotoc {

///
/// @class ImpactStage
/// @brief Impact stage computations for optimal control problems.
///
class ImpactStage {
public:
  ///
  /// @brief Constructs the impact stage.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  ///
  ImpactStage(const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints,
              const std::shared_ptr<ContactSequence>& contact_sequence);

  ///
  /// @brief Default constructor.  
  ///
  ImpactStage();

  ///
  /// @brief Default destructor. 
  ///
  ~ImpactStage() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ImpactStage(const ImpactStage&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpactStage& operator=(const ImpactStage&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpactStage(ImpactStage&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpactStage& operator=(ImpactStage&&) noexcept = default;

  ///
  /// @brief Creates the data.
  /// @param[in] robot Robot model. 
  ///
  OCPData createData(const Robot& robot) const;

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  bool isFeasible(Robot& robot, const GridInfo& grid_info,
                  const SplitSolution& s, OCPData& data) const;

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  void initConstraints(Robot& robot, const GridInfo& grid_info, 
                       const SplitSolution& s, OCPData& data) const;

  ///
  /// @brief Computes the stage cost and constraint violation.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] kkt_residual KKT residual of this stage. 
  ///
  void evalOCP(Robot& robot, const GridInfo& grid_info, const SplitSolution& s, 
               const SplitSolution& s_next, OCPData& data,
               SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the KKT residual and matrix of this stage.
  /// @param[in, out] robot Robot model. 
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in] q_prev Configuration at the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] kkt_matrix KKT matrix of this stage. 
  /// @param[in, out] kkt_residual KKT residual of this stage. 
  ///
  void evalKKT(Robot& robot, const GridInfo& grid_info, 
               const Eigen::VectorXd& q_prev, const SplitSolution& s, 
               const SplitSolution& s_next, OCPData& data,
               SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Expands the condensed primal variables, i.e., computes the Newton 
  /// direction of the condensed primal variables of this stage.
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in, out] d Split direction of this stage.
  /// 
  void expandPrimal(const GridInfo& grid_info, OCPData& data, SplitDirection& d) const;

  ///
  /// @brief Expands the condensed dual variables, i.e., computes the Newton 
  /// direction of the condensed dual variables of this stage.
  /// @param[in] grid_info Grid info of this stage.
  /// @param[in, out] data Data of this stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void expandDual(const GridInfo& grid_info, OCPData& data, 
                  const SplitDirection& d_next, 
                  SplitDirection& d) const;

  ///
  /// @brief Returns maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  /// @param[in] data Data of this stage. 
  /// @return Maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  ///
  double maxPrimalStepSize(const OCPData& data) const;

  ///
  /// @brief Returns maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  /// @param[in] data Data of this stage. 
  /// @return Maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  ///
  double maxDualStepSize(const OCPData& data) const;

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  /// @param[in, out] data Data of this stage. 
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s, OCPData& data) const;

  ///
  /// @brief Updates dual variables of this stage.
  /// @param[in] dual_step_size Dual step size. 
  /// @param[in, out] data Data of this stage. 
  ///
  void updateDual(const double dual_step_size, OCPData& data) const;

private:
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  std::shared_ptr<ContactSequence> contact_sequence_;
};

} // namespace robotoc

#endif // ROBOTOC_IMPACT_STAGE_HPP_