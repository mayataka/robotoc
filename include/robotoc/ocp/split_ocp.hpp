#ifndef ROBOTOC_SPLIT_OCP_HPP_
#define ROBOTOC_SPLIT_OCP_HPP_

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
#include "robotoc/ocp/state_equation.hpp"
#include "robotoc/ocp/contact_dynamics.hpp"
#include "robotoc/ocp/switching_constraint.hpp"
#include "robotoc/ocp/switching_constraint_residual.hpp"
#include "robotoc/ocp/switching_constraint_jacobian.hpp"


namespace robotoc {

///
/// @class SplitOCP
/// @brief An optimal control problem for Riccati recursion algorithm split
/// into a time stage. 
///
class SplitOCP {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
           const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  SplitOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitOCP(const SplitOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitOCP& operator=(const SplitOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitOCP(SplitOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitOCP& operator=(SplitOCP&&) noexcept = default;

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
  /// @brief Initializes the constraints, i.e., copies the slack and dual 
  /// variables from another time stage. 
  /// @param[in] other Split optimal control problem at another time stage.
  ///
  void initConstraints(const SplitOCP& other);

  ///
  /// @brief Gets the const reference to the constraints data. 
  /// @return const reference to the constraints data. 
  ///
  const ConstraintsData& constraintsData() const;

  ///
  /// @brief Computes the stage cost and constraint violation.
  /// Used in the line search.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in] q_next Configuration at the next time stage.
  /// @param[in] v_next Generaized velocity at the next time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void evalOCP(Robot& robot, const ContactStatus& contact_status,
               const double t, const double dt, const SplitSolution& s, 
               const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next,
               SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the stage cost and constraint violation.
  /// Used in the line search.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in] q_next Configuration at the next time stage.
  /// @param[in] v_next Generaized velocity at the next time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  /// @param[in] dt_next Time step of the next time stage. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  ///
  void evalOCP(Robot& robot, const ContactStatus& contact_status,
               const double t, const double dt, const SplitSolution& s, 
               const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next,
               SplitKKTResidual& kkt_residual,
               const ImpulseStatus& impulse_status, const double dt_next, 
               SwitchingConstraintResidual& sc_residual);

  ///
  /// @brief Computes the KKT residual of this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  template <typename SplitSolutionType>
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status,
                          const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolutionType& s_next,
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual of this time stage including the 
  /// switching constraint.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  /// @param[in] dt_next Time step of the next time stage. 
  /// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  ///
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status,
                          const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolution& s_next, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual, 
                          const ImpulseStatus& impulse_status, 
                          const double dt_next, 
                          SwitchingConstraintJacobian& sc_jacobian,
                          SwitchingConstraintResidual& sc_residual);

  ///
  /// @brief Computes the KKT system of this time stage, i.e., the condensed
  /// KKT matrix and KKT residual of this time stage for Newton's method.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  template <typename SplitSolutionType>
  void computeKKTSystem(Robot& robot, const ContactStatus& contact_status, 
                        const double t, const double dt, 
                        const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                        const SplitSolutionType& s_next, 
                        SplitKKTMatrix& kkt_matrix,
                        SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT system of this time stage, i.e., the condensed
  /// KKT matrix and KKT residual of this time stage including the 
  /// switching constraint for Newton's method.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  /// @param[in] dt_next Time step of the next time stage. 
  /// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  ///
  void computeKKTSystem(Robot& robot, const ContactStatus& contact_status, 
                        const double t, const double dt, 
                        const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                        const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix,
                        SplitKKTResidual& kkt_residual, 
                        const ImpulseStatus& impulse_status, 
                        const double dt_next, 
                        SwitchingConstraintJacobian& sc_jacobian,
                        SwitchingConstraintResidual& sc_residual);

  ///
  /// @brief Computes the initial state direction using the result of  
  /// SplitOCP::computeKKTSystem().
  /// @param[in] robot Robot model. 
  /// @param[in] q0 Initial configuration. 
  /// @param[in] v0 Initial generalized velocity. 
  /// @param[in] s0 Split solution at the initial stage. 
  /// @param[in] d0 Split direction at the initial stage. 
  ///
  void computeInitialStateDirection(const Robot& robot, 
                                    const Eigen::VectorXd& q0, 
                                    const Eigen::VectorXd& v0, 
                                    const SplitSolution& s0, 
                                    SplitDirection& d0) const;

  ///
  /// @brief Expands the condensed primal variables, i.e., computes the Newton 
  /// direction of the condensed primal variables of this stage.
  /// @param[in] dt Time step of this time stage. This is used only when sto is 
  /// true.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// @param[in] sto If true, the sensitivity w.r.t. the switching time is 
  /// considered. If false, it is not considered. 
  /// 
  void expandPrimal(const double dt, const SplitSolution& s, 
                    SplitDirection& d, const bool sto);

  ///
  /// @brief Expands the condensed dual variables, i.e., computes the Newton 
  /// direction of the condensed dual variables of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] d_next Split direction of the next time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// @param[in] sto If true, the sensitivity w.r.t. the switching time is 
  /// considered. If false, it is not considered. 
  /// 
  template <typename SplitDirectionType>
  void expandDual(const double dt, const SplitDirectionType& d_next, 
                  SplitDirection& d, const bool sto);

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
  /// function, SplitOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of this time stage.
  /// @param[in] dt Time step of this time stage.
  /// @return The squared norm of the kKT residual.
  ///
  double KKTError(const SplitKKTResidual& kkt_residual, const double dt) const;

  ///
  /// @brief Returns the stage cost of this time stage for the line search.
  /// Before calling this function, SplitOCP::evalOCP(), 
  /// SplitOCP::computeKKTResidual(), SplitOCP::computeKKTSystem() must be called.
  /// @return Stage cost of this time stage.
  /// 
  double stageCost() const;

  ///
  /// @brief Returns the constraint violation of this time stage for the 
  /// line search. 
  /// Before calling this function, SplitOCP::evalOCP(), 
  /// SplitOCP::computeKKTResidual(), SplitOCP::computeKKTSystem() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @param[in] dt Time step of this time stage. 
  /// @return The constraint violation of this time stage.
  ///
  double constraintViolation(const SplitKKTResidual& kkt_residual, 
                             const double dt) const;

  ///
  /// @brief Returns the constraint violation of this time stage for the 
  /// line search. 
  /// Before calling this function, SplitOCP::evalOCP(), 
  /// SplitOCP::computeKKTResidual(), SplitOCP::computeKKTSystem() must be called.
  /// @param[in] kkt_residual KKT residual of this impulse stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] sc_residual Residual of the switching constraint. 
  /// @return The constraint violation of this time stage.
  ///
  double constraintViolation(
      const SplitKKTResidual& kkt_residual, const double dt,
      const SwitchingConstraintResidual& sc_residual) const;

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  StateEquation state_equation_;
  ContactDynamics contact_dynamics_;
  double stage_cost_;

};

} // namespace robotoc

#include "robotoc/ocp/split_ocp.hxx"

#endif // ROBOTOC_SPLIT_OCP_HPP_