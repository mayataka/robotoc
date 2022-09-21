#ifndef ROBOTOC_TERMINAL_OCP_HPP_
#define ROBOTOC_TERMINAL_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

///
/// @class TerminalOCP 
/// @brief An optimal control problem for Riccati recursion algorithm split
/// into a terminal stage. 
///
class TerminalOCP {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  TerminalOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  TerminalOCP();
  
  ///
  /// @brief Destructor. 
  ///
  ~TerminalOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  TerminalOCP(const TerminalOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TerminalOCP& operator=(const TerminalOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TerminalOCP(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TerminalOCP& operator=(TerminalOCP&&) noexcept = default;

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of the terminal stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. 
  /// @param[in] time_stage Time stage.
  /// @param[in] s Split solution of the terminal stage.
  ///
  void initConstraints(Robot& robot, const int time_stage, 
                       const SplitSolution& s);

  ///
  /// @brief Computes the terminal cost and constraint violation.
  /// Used in the line search.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info of this time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void evalOCP(Robot& robot, const GridInfo& grid_info, const SplitSolution& s, 
               SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual of the terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info of this time stage.
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of the terminal stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of the terminal stage.
  /// @param[in, out] kkt_residual Split KKT residual of the terminal stage.
  ///
  void computeKKTResidual(Robot& robot, const GridInfo& grid_info,
                          const Eigen::VectorXd& q_prev, const SplitSolution& s,
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT system of the terminal stage, i.e., the condensed
  /// KKT matrix and KKT residual of the terminal stage for Newton's method.
  /// @param[in] robot Robot model. 
  /// @param[in] grid_info Grid info of this time stage.
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of the terminal stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of the terminal stage.
  /// @param[in, out] kkt_residual Split KKT residual of the terminal stage.
  ///
  void computeKKTSystem(Robot& robot, const GridInfo& grid_info,
                        const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                        SplitKKTMatrix& kkt_matrix, 
                        SplitKKTResidual& kkt_residual);

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
  /// @brief Expands the condensed primal variables, i.e., computes the Newton 
  /// direction of the condensed primal variables of the terminal stage.
  /// @param[in] s Split solution of the terminal stage.
  /// @param[in, out] d Split direction of the terminal stage.
  /// 
  void expandPrimal(const SplitSolution& s, SplitDirection& d);

  ///
  /// @brief Expands the condensed dual variables, i.e., computes the Newton 
  /// direction of the condensed dual variables of the terminal stage.
  /// @param[in, out] d Split direction of the terminal stage.
  /// 
  void expandDual(SplitDirection& d);

  ///
  /// @brief Updates primal variables of the terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] d Split direction of the terminal stabe.
  /// @param[in, out] s Split solution of the terminal stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s) const;

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size.
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Returns the KKT residual of the terminal stage. Before calling this 
  /// function, TerminalOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of the terminal stage.
  /// @return The squared norm of the kKT residual.
  ///
  static double KKTError(const SplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the terminal cost for the line search.
  /// Before calling this function, TerminalOCP::evalOCP(), 
  /// TerminalOCP::computeKKTResidual(), or TerminalOCP::computeKKTSystem() 
  /// must be called.
  /// @param[in] include_cost_barrier If true, includes the cost due to the 
  /// barrier function. Default is true.
  /// 
  double terminalCost(const bool include_cost_barrier=true) const;

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  TerminalStateEquation state_equation_;
  double terminal_cost_, barrier_cost_;

};

} // namespace robotoc

#endif // IDOCPT_TERMINAL_OCP_HPP_