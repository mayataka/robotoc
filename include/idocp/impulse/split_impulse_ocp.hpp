#ifndef IDOCP_SPLIT_IMPULSE_OCP_HPP_
#define IDOCP_SPLIT_IMPULSE_OCP_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"
#include "idocp/impulse/impulse_riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

///
/// @class SplitImpulseOCP
/// @brief Split optimal control problem of a single stage. 
///
class SplitImpulseOCP {
public:
  ///
  /// @brief Construct a split optimal control problem.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitImpulseOCP(const Robot& robot, 
                  const std::shared_ptr<ImpulseCostFunction>& cost,
                  const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor.  
  ///
  SplitImpulseOCP();

  ///
  /// @brief Destructor. 
  ///
  ~SplitImpulseOCP();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitImpulseOCP(const SplitImpulseOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitImpulseOCP& operator=(const SplitImpulseOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitImpulseOCP(SplitImpulseOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitImpulseOCP& operator=(SplitImpulseOCP&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  void initConstraints(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Linearize the OCP for Newton's method around the current solution, 
  /// i.e., computes the KKT residual and Hessian.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void linearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                    const double t, const Eigen::VectorXd& q_prev, 
                    const ImpulseSplitSolution& s, 
                    const SplitSolution& s_next);

  ///
  /// @brief Computes the Riccati factorization of this stage from the 
  /// factorization of the previous stage.
  /// @param[in] riccati_next Riccati factorization of the next stage.
  /// @param[out] riccati Riccati factorization of this stage.
  /// 
  void backwardRiccatiRecursion(const RiccatiFactorization& riccati_next,
                                RiccatiFactorization& riccati);

  ///
  /// @brief Computes the Newton direction of the state of this stage from the 
  /// one of the previous stage.
  /// @param[in] d Split direction of this stage.
  /// @param[in] d_next Split direction of the next stage.
  /// 
  void forwardRiccatiRecursion(ImpulseSplitDirection& d,
                               SplitDirection& d_next);

  ///
  /// @brief Computes the Newton direction of the condensed variables of this 
  /// stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedDirection(Robot& robot, const ImpulseSplitSolution& s, 
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
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage. The stage cost is recomputed. The violation of the  
  /// constriants is not computed. Instead, the previously computed residual  
  /// computed by SplitOCP::linearizeOCP or 
  /// SplitOCP::computeKKTResidual, is used.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndConstraintViolation(
      Robot& robot, const double t, const ImpulseSplitSolution& s);

  ///Impulse
  /// @brief Returns the stage cost and L1-norm of the violation of constraints 
  /// of this stage under step_size. The split solution of this stage and the 
  /// state of the next stage are computed by step_size temporary. 
  /// The stage cost and the violation of the constriants are computed based on
  /// the temporary solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] step_size Step size for the primal variables. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[in] d_next Split direction of the next stage.
  /// @return The stage cost and L1-norm of the constraints violation.
  ///
  std::pair<double, double> costAndConstraintViolation(
      Robot& robot, const ContactStatus& contact_status, const double step_size, 
      const double t, const ImpulseSplitSolution& s, 
      const ImpulseSplitDirection& d, 
      const SplitSolution& s_next, const SplitDirection& d_next);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] step_size Dula step size of the OCP. 
  ///
  void updateDual(const double step_size);

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] step_size Primal step size of the OCP. 
  /// @param[in] riccati Riccati factorization of this stage.
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double step_size, 
                    const RiccatiFactorization& riccati, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void computeKKTResidual(Robot& robot, const ContactStatus& contact_status,
                          const double t, const Eigen::VectorXd& q_prev, 
                          const ImpulseSplitSolution& s, 
                          const SplitSolution& s_next);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<ImpulseCostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseKKTResidual kkt_residual_;
  ImpulseKKTMatrix kkt_matrix_;
  ImpulseDynamicsForwardEuler impulse_dynamics_;
  ImpulseRiccatiMatrixFactorizer riccati_factorizer_;
  ImpulseSplitSolution s_tmp_; /// @brief Temporary split solution used in line search.
  int dimv_, dimf_, dimc_;
  double stage_cost_, constraint_violation_;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  inline void setContactStatusForKKT(const ContactStatus& contact_status) {
    kkt_residual_.setContactStatus(contact_status);
    kkt_matrix_.setContactStatus(contact_status);
    dimf_ = contact_status.dimf();
    dimc_ = contact_status.dimf();
  }

  double cost(Robot& robot, const double t, const ImpulseSplitSolution& s);

  double constraintViolation() const;

};

} // namespace idocp


#endif // IDOCP_SPLIT_IMPULSE_OCP_HPP_ 