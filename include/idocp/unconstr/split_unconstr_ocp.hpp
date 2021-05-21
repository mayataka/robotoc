#ifndef IDOCP_SPLIT_UNCONSTR_OCP_HPP_
#define IDOCP_SPLIT_UNCONSTR_OCP_HPP_

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
#include "idocp/unconstr/unconstr_dynamics.hpp"


namespace idocp {

///
/// @class SplitUnconstrOCP
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// Riccati recursion algorithm split into a time stage. 
///
class SplitUnconstrOCP {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitUnconstrOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  SplitUnconstrOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnconstrOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnconstrOCP(const SplitUnconstrOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitUnconstrOCP& operator=(const SplitUnconstrOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitUnconstrOCP(SplitUnconstrOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnconstrOCP& operator=(SplitUnconstrOCP&&) noexcept = default;

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
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix,
                    SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed variables  
  /// at this stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedDirection(Robot& robot, const double dt, 
                                 const SplitSolution& s, 
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
  /// @brief Computes the KKT residual of this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void computeKKTResidual(Robot& robot, const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolution& s_next, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of this time stage. Before calling this 
  /// function, SplitUnconstrOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in] dt Time step of this time stage.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const SplitKKTResidual& kkt_residual, 
                                const double dt) const;

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
  /// @param[in] s Split solution of this time stage.
  /// @param[in] q_next Configuration at the next time stage.
  /// @param[in] v_next Generaized velocity at the next time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @return Constraint violation of this time stage.
  ///
  double constraintViolation(Robot& robot,  const double t, const double dt, 
                             const SplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next,
                             SplitKKTResidual& kkt_residual);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  UnconstrDynamics unconstr_dynamics_;
  bool use_kinematics_;
  double stage_cost_, constraint_violation_;

};

} // namespace idocp

#include "idocp/unconstr/split_unconstr_ocp.hxx"

#endif // IDOCP_SPLIT_UNCONSTR_OCP_HPP_ 