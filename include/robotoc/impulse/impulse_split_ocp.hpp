#ifndef ROBOTOC_IMPULSE_SPLIT_OCP_HPP_
#define ROBOTOC_IMPULSE_SPLIT_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/cost_function_data.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/impulse/impulse_state_equation.hpp"
#include "robotoc/impulse/impulse_dynamics.hpp"
#include "robotoc/ocp/split_direction.hpp"


namespace robotoc {

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
  /// @brief Computes the impulse stage cost and constraint violation.
  /// Used in the line search.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] q_next Configuration at the next time stage.
  /// @param[in] v_next Generaized velocity at the next time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void evalOCP(Robot& robot, const ImpulseStatus& impulse_status,
               const double t, const ImpulseSplitSolution& s, 
               const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next,
               ImpulseSplitKKTResidual& kkt_residual);

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
  /// @brief Computes the KKT system of this impuse stage, i.e., the condensed
  /// KKT matrix and KKT residual of this impulse stage for Newton's method.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status of this impulse stage. 
  /// @param[in] t Time of this impulse stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage.
  ///
  void computeKKTSystem(Robot& robot, const ImpulseStatus& impulse_status, 
                        const double t, const Eigen::VectorXd& q_prev, 
                        const ImpulseSplitSolution& s, 
                        const SplitSolution& s_next, 
                        ImpulseSplitKKTMatrix& kkt_matrix, 
                        ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Expands the condensed primal variables, i.e., computes the Newton 
  /// direction of the condensed primal variables of this impulse stage.
  /// @param[in] s Split solution of this impulse stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void expandPrimal(const ImpulseSplitSolution& s, ImpulseSplitDirection& d);

  ///
  /// @brief Expands the condensed dual variables, i.e., computes the Newton 
  /// direction of the condensed dual variables of this impulse stage.
  /// @param[in] d_next Split direction of the next time stage.
  /// @param[in, out] d Split direction of this impulse stage.
  /// 
  void expandDual(const SplitDirection& d_next, ImpulseSplitDirection& d);

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
  /// @brief Returns the KKT residual of this impulse stage. Before calling this 
  /// function, ImpulseSplitOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @return The squared norm of the kKT residual.
  ///
  double KKTError(const ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Returns the stage cost of this impulse stage for the line search.
  /// Before calling this function, ImpulseSplitOCP::evalOCP(),
  /// ImpulseSplitOCP::computeKKTResidual(),
  /// or ImpulseSplitOCP::computeKKTSystem() must be called.
  /// @return The stage cost of this impulse stage.
  /// 
  double stageCost() const;

  ///
  /// @brief Returns the constraint violation of this impulse stage for the 
  /// line search.
  /// Before calling this function, ImpulseSplitOCP::evalOCP(),
  /// ImpulseSplitOCP::computeKKTResidual(),
  /// or ImpulseSplitOCP::computeKKTSystem() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @return The constraint violation of this impulse stage.
  ///
  double constraintViolation(const ImpulseSplitKKTResidual& kkt_residual) const;

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseStateEquation state_equation_;
  ImpulseDynamics impulse_dynamics_;
  double stage_cost_;

};

} // namespace robotoc

#include "robotoc/impulse/impulse_split_ocp.hxx"

#endif // ROBOTOC_IMPULSE_SPLIT_OCP_HPP_ 