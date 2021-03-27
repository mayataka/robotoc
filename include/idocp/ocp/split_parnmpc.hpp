#ifndef IDOCP_SPLIT_PARNMPC_HPP_ 
#define IDOCP_SPLIT_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

///
/// @class SplitParNMPC
/// @brief An optimal control problem for ParNMPC algorithm split into a time 
/// stage. 
///
class SplitParNMPC {
public:
  ///
  /// @brief Constructs a split optimal control problem.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
               const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  SplitParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitParNMPC(const SplitParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitParNMPC& operator=(const SplitParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitParNMPC(SplitParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitParNMPC& operator=(SplitParNMPC&&) noexcept = default;

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
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  template <typename SplitSolutionType>
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
                    const SplitSolution& s, const SplitSolutionType& s_next, 
                    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian. Also linearize the 
  /// switching constraints.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  ///
  template <typename SplitSolutionType>
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const double t, const double dt, 
                    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
                    const SplitSolution& s, const SplitSolutionType& s_next, 
                    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual,
                    const ImpulseStatus& impulse_status);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot, const double dt, 
                                       const SplitSolution& s, 
                                       SplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this time stage.
  /// @param[in] robot Robot model. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] kkt_matrix KKT matrix of this time stage.
  /// @param[in, out] kkt_residual KKT residual of this time stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void computeCondensedDualDirection(const Robot& robot, const double dt, 
                                     const SplitKKTMatrix& kkt_matrix, 
                                     SplitKKTResidual& kkt_residual,
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
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  template <typename SplitSolutionType>
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                          const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s, 
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
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in] s_next Split solution of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  ///
  template <typename SplitSolutionType>
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                          const double t, const double dt, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                          const SplitSolutionType& s_next, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual,
                          const ImpulseStatus& impulse_status);

  ///
  /// @brief Returns the KKT residual of this time stage. Before calling this 
  /// function, SplitOCP::computeKKTResidual() must be called.
  /// @param[in] kkt_residual KKT residual of this time stage.
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
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual KKT residual of this time stage.
  /// @return Constraint violation of this time stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dt, 
                             const Eigen::VectorXd& q_prev,
                             const Eigen::VectorXd& v_prev,
                             const SplitSolution& s,
                             SplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the constraint violation of this time stage including
  /// the switching constraint for line search.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] t Time of this time stage. 
  /// @param[in] dt Time of this time stage. 
  /// @param[in] q_prev Configuration at the previous time stage.
  /// @param[in] v_prev Generalized velocity at the previous time stage.
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual KKT residual of this time stage.
  /// @param[in] impulse_status Impulse status at the switching instant. 
  /// @return Constraint violation of this time stage.
  ///
  double constraintViolation(Robot& robot, const ContactStatus& contact_status, 
                             const double t, const double dt, 
                             const Eigen::VectorXd& q_prev,
                             const Eigen::VectorXd& v_prev,
                             const SplitSolution& s,
                             SplitKKTResidual& kkt_residual,
                             const ImpulseStatus& impulse_status);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ContactDynamics contact_dynamics_;
  bool use_kinematics_, has_floating_base_;
  double stage_cost_, constraint_violation_;

};

} // namespace idocp

#include "idocp/ocp/split_parnmpc.hxx"

#endif // IDOCP_SPLIT_PARNMPC_HPP_