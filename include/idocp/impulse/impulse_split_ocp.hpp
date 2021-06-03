#ifndef IDOCP_IMPULSE_SPLIT_OCP_HPP_
#define IDOCP_IMPULSE_SPLIT_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"
#include "idocp/impulse/impulse_dynamics.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class ImpulseSplitOCP
/// @brief An optimal control problem for Riccati recursion algorithm split
/// into an impulse stage. 
///
class ImpulseSplitOCP {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  ImpulseSplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                  const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  ImpulseSplitOCP();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitOCP(const ImpulseSplitOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseSplitOCP& operator=(const ImpulseSplitOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitOCP(ImpulseSplitOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitOCP& operator=(ImpulseSplitOCP&&) noexcept = default;

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this impulse stage.
  ///
  bool isFeasible(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Initializes the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. 
  /// @param[in] s Split solution of this impulse stage.
  ///
  void initConstraints(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Linearizes the split optimal control problem for Newton's method 
  /// around the current solution, i.e., computes the KKT residual and Hessian.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void linearizeOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                    const double t, const Eigen::VectorXd& q_prev, 
                    const ImpulseSplitSolution& s, 
                    const SplitSolution& s_next, 
                    ImpulseSplitKKTMatrix& kkt_matrix, 
                    ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this impulse stage.
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void computeCondensedPrimalDirection(const ImpulseSplitSolution& s, 
                                       ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this impulse stage.
  /// @param[in] kkt_matrix KKT matrix of this impulse stage.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @param[in] d_next Split direction of the next time stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void computeCondensedDualDirection(const ImpulseSplitKKTMatrix& kkt_matrix, 
                                     const ImpulseSplitKKTResidual& kkt_residual,
                                     const SplitDirection& d_next, 
                                     ImpulseSplitDirection& d);

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
  /// @brief Updates primal variables of this impulse stage.
  /// @param[in] robot Robot model. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] d Split direction of this impulse stage.
  /// @param[in, out] s Split solution of this impulse stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size.
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of this impulse stage.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void computeKKTResidual(Robot& robot, const ImpulseStatus& impulse_status,
                          const double t, const Eigen::VectorXd& q_prev, 
                          const ImpulseSplitSolution& s, 
                          const SplitSolution& s_next,
                          ImpulseSplitKKTMatrix& kkt_matrix, 
                          ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of this impulse stage. Before calling this 
  /// function, ImpulseSplitOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the stage cost of this impulse stage for line search.
  /// @param[in] robot Robot model. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] primal_step_size Primal step size. Default is 0.
  /// @return Stage cost of this impulse stage.
  /// 
  double stageCost(Robot& robot, const double t, const ImpulseSplitSolution& s, 
                   const double primal_step_size=0);

  ///
  /// @brief Computes the constraint violation of this impulse stage for line 
  /// search.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status oif this impulse stage. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] q_next Configuration at the next time stage.
  /// @param[in] v_next Generaized velocity at the next time stage.
  /// @param[in, out] kkt_residual KKT residual of this impulse stage.
  /// @return Constraint violation of this impulse stage.
  ///
  double constraintViolation(Robot& robot, const ImpulseStatus& impulse_status, 
                             const double t, const ImpulseSplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next,
                             ImpulseSplitKKTResidual& kkt_residual);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseStateEquation state_equation_;
  ImpulseDynamics impulse_dynamics_;
  double impulse_cost_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_ocp.hxx"

#endif // IDOCP_IMPULSE_SPLIT_OCP_HPP_ 