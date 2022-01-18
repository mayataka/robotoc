#ifndef ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_
#define ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/unconstr/unconstr_state_equation.hpp"
#include "robotoc/unconstr/unconstr_dynamics.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

///
/// @class SplitUnconstrParNMPC
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// ParNMPC algorithm split into a time stage. 
///
class SplitUnconstrParNMPC {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitUnconstrParNMPC(const Robot& robot, 
                       const std::shared_ptr<CostFunction>& cost,
                       const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  SplitUnconstrParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnconstrParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnconstrParNMPC(const SplitUnconstrParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitUnconstrParNMPC& operator=(const SplitUnconstrParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitUnconstrParNMPC(SplitUnconstrParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnconstrParNMPC& operator=(SplitUnconstrParNMPC&&) noexcept = default;

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this time stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. 
  /// @param[in] time_stage Time stage.
  /// @param[in] s Split solution of this time stage.
  ///
  void initConstraints(Robot& robot, const int time_stage, 
                       const SplitSolution& s);

  ///
  /// @brief Computes the stage cost and constraint violation.
  /// Used in the line search.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void evalOCP(Robot& robot, const GridInfo& grid_info, 
               const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
               const SplitSolution& s, SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual of this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void computeKKTResidual(Robot& robot, const GridInfo& grid_info,
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s,
                          const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT system of this time stage, i.e., the condensed
  /// KKT matrix and KKT residual of this time stage for Newton's method.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void computeKKTSystem(Robot& robot, const GridInfo& grid_info,
                        const Eigen::VectorXd& q_prev, 
                        const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                        const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
                        SplitKKTResidual& kkt_residual);

  ///
  /// @brief Expands the primal and dual variables, i.e., computes the Newton 
  /// direction of the condensed variables of this stage.
  /// @param[in] dt Time step of this time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void expandPrimalAndDual(const double dt, const SplitSolution& s, 
                           const SplitKKTMatrix& kkt_matrix,
                           const SplitKKTResidual& kkt_residual,
                           SplitDirection& d);

  ///
  /// @brief Returns maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the primal variables that satisfies 
  /// the inequality constraints.
  ///
  double maxPrimalStepSize();

  ///
  /// @brief Returns maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  /// @return Maximum stap size of the dual variables that satisfies 
  /// the inequality constraints.
  ///
  double maxDualStepSize();

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Returns the KKT residual of this time stage. Before calling this 
  /// function, computeKKTResidual() must be called.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in] dt Time step of this time stage.
  /// @return The squared norm of the kKT residual.
  ///
  double KKTError(const SplitKKTResidual& kkt_residual, const double dt) const;

  ///
  /// @brief Returns the stage cost of this time stage for the line search.
  /// Before calling this function, SplitUnconstrParNMPC::evalOCP(), 
  /// SplitUnconstrParNMPC::computeKKTResidual(), or
  /// SplitUnconstrParNMPC::computeKKTSystem() must be called.
  /// @return Stage cost of this time stage.
  /// 
  double stageCost() const;

  ///
  /// @brief Returns the constraint violation of this time stage for the 
  /// line search. Before calling this function, 
  /// SplitUnconstrParNMPC::evalOCP(), 
  /// SplitUnconstrParNMPC::computeKKTResidual(), or
  /// SplitUnconstrParNMPC::computeKKTSystem() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @param[in] dt Time step of this time stage. 
  /// @return The constraint violation of this time stage.
  ///
  double constraintViolation(const SplitKKTResidual& kkt_residual, 
                             const double dt) const;

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  UnconstrDynamics unconstr_dynamics_;
  ContactStatus contact_status_;
  bool use_kinematics_;
  double stage_cost_;

};

} // namespace robotoc

#endif // ROBOTOC_SPLIT_UNCONSTR_PARNMPC_HPP_ 