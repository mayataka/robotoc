#ifndef IDOCP_TERMINAL_UNPARNMPC_HPP_
#define IDOCP_TERMINAL_UNPARNMPC_HPP_

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
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/unconstrained_dynamics.hpp"


namespace idocp {

///
/// @class TerminalParNMPC
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// ParNMPC algorithm split into the terminal stage. 
///
class TerminalUnParNMPC {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  TerminalUnParNMPC(const Robot& robot, 
                    const std::shared_ptr<CostFunction>& cost,
                    const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  TerminalUnParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~TerminalUnParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  TerminalUnParNMPC(const TerminalUnParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  TerminalUnParNMPC& operator=(const TerminalUnParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  TerminalUnParNMPC(TerminalUnParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  TerminalUnParNMPC& operator=(TerminalUnParNMPC&&) noexcept = default;

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
  void initConstraints(Robot& robot, const int time_step, 
                       const SplitSolution& s);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] unkkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] unkkt_residual Split KKT residual of this time stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
                    const SplitSolution& s, SplitUnKKTMatrix& unkkt_matrix,
                    SplitUnKKTResidual& unkkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed variables of 
  /// this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedDirection(Robot& robot, const double dt, 
                                 const SplitSolution& s, SplitDirection& d);

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
  /// @brief Computes the KKT residual of this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  ///
  void computeKKTResidual(Robot& robot, const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s);

  ///
  /// @brief Returns the KKT residual of this time stage. Before calling this 
  /// function, SplitOCP::computeKKTResidual() must be called.
  /// @param[in] dt Time step of this time stage.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const double dt) const;

  ///
  /// @brief Computes the stage cost of this time stage for line search.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in] primal_step_size Primal step size. Default is 0.
  /// @return Stage cost of this time stage.
  /// 
  double stageCost(Robot& robot, const double t, const double dt, 
                   const SplitSolution& s, const double primal_step_size=0);

  ///
  /// @brief Computes the constraint violation of this time stage for line 
  /// search.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @return Constraint violation of this time stage.
  ///
  double constraintViolation(Robot& robot,  const double t, const double dt, 
                             const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const SplitSolution& s);

  ///
  /// @brief Computes the terminal cost Hessian to initialize the auxiliary 
  /// matrices for backward correction. 
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] Qxx Hessian with respect to the state.
  ///
  template <typename MatrixType>
  void computeTerminalCostHessian(Robot& robot, const double t, 
                                  const SplitSolution& s, 
                                  const Eigen::MatrixBase<MatrixType>& Qxx);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  UnconstrainedDynamics unconstrained_dynamics_;
  bool use_kinematics_;
  SplitKKTMatrix kkt_matrix_;
  SplitKKTResidual kkt_residual_;
  double stage_cost_, constraint_violation_;

};

} // namespace idocp

#include "idocp/unocp/terminal_unparnmpc.hxx"

#endif // IDOCP_TERMINAL_UNPARNMPC_HPP_ 