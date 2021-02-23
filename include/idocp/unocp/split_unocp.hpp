#ifndef IDOCP_SPLIT_UNOCP_HPP_
#define IDOCP_SPLIT_UNOCP_HPP_

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
/// @class SplitUnOCP
/// @brief Split unconstrained optimal control problem of a single stage. 
///
class SplitUnOCP {
public:
  ///
  /// @brief Construct a split optimal control problem.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitUnOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
             const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  SplitUnOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnOCP(const SplitUnOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitUnOCP& operator=(const SplitUnOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitUnOCP(SplitUnOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnOCP& operator=(SplitUnOCP&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const SplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const int time_step, 
                       const SplitSolution& s);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] unkkt_matrix Condensed KKT matrix of this stage.
  /// @param[out] unkkt_residual Condensed KKT residual of this stage.
  ///
  void linearizeOCP(Robot& robot, const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                    const SplitSolution& s_next, SplitUnKKTMatrix& unkkt_matrix,
                    SplitUnKKTResidual& unkkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed variables  
  /// at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] d Split direction of this stage.
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
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] primal_step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const SplitDirection& d, SplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size of the OCP. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @tparam SplitSolutionType Type of the split solution at the next stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void computeKKTResidual(Robot& robot, const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolution& s_next);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitUnOCP::linearizeOCP() or 
  /// SplitUnOCP::computeKKTResidual() must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const double dt) const;

  ///
  /// @brief Computes the stage cost of this stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] primal_step_size Primal step size of the OCP. Default is 0.
  /// @return Stage cost of this stage.
  /// 
  double stageCost(Robot& robot, const double t, const double dt, 
                   const SplitSolution& s, const double primal_step_size=0);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] dt Length of the discretization of the horizon.
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_next Configuration at the next stage.
  /// @param[in] v_next Generaized velocity at the next stage.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot,  const double t, const double dt, 
                             const SplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

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

#include "idocp/unocp/split_unocp.hxx"

#endif // IDOCP_SPLIT_UNOCP_HPP_