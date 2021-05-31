#ifndef IDOCP_TERMINAL_OCP_HPP_
#define IDOCP_TERMINAL_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"


namespace idocp {

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
  /// @param[in] s Split solution of this terminal stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. 
  /// @param[in] time_stage Time stage.
  /// @param[in] s Split solution of this terminal stage.
  ///
  void initConstraints(Robot& robot, const int time_stage, 
                       const SplitSolution& s);

  ///
  /// @brief Linearizes the split optimal control problem for Newton's method 
  /// around the current solution, i.e., computes the KKT residual and Hessian.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this terminal stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this terminal stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this terminal stage.
  /// @param[in, out] kkt_residual Split KKT residual of this terminal stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const Eigen::VectorXd& q_prev, 
                    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
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
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this terminal stage.
  /// @param[in, out] d Split direction of this terminal stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot, const SplitSolution& s, 
                                       SplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] kkt_matrix KKT matrix of this terminal stage.
  /// @param[in, out] kkt_residual KKT residual of this terminal stage.
  /// @param[in, out] d Split direction of this terminal stage.
  /// 
  void computeCondensedDualDirection(const Robot& robot, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual,
                                     SplitDirection& d);

  ///
  /// @brief Updates primal variables of this terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] d Split direction of this terminal stabe.
  /// @param[in, out] s Split solution of this terminal stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s) const;

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size.
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of this terminal stage.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this terminal stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this terminal stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this terminal stage.
  /// @param[in, out] kkt_residual Split KKT residual of this terminal stage.
  ///
  void computeKKTResidual(Robot& robot, const double t, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s,
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of this terminal stage. Before calling this 
  /// function, TerminalOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of this terminal stage.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the terminal cost of this terminal stage for line search.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this terminal stage. 
  /// @param[in] s Split solution of this terminal stage.
  /// @return Terminal cost of this terminal stage.
  /// 
  double terminalCost(Robot& robot, const double t, const SplitSolution& s);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  bool use_kinematics_;
  double terminal_cost_;

};

} // namespace idocp

#include "idocp/ocp/terminal_ocp.hxx"

#endif // IDOCPT_TERMINAL_OCP_HPP_