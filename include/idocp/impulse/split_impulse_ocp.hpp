#ifndef IDOCP_SPLIT_IMPULSE_OCP_HPP_
#define IDOCP_SPLIT_IMPULSE_OCP_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"


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
  /// @param[in] cost Shared ptr to the impulse cost function.
  /// @param[in] constraints Shared ptr to the impulse constraints.
  ///
  SplitImpulseOCP(const Robot& robot, 
                  const std::shared_ptr<ImpulseCostFunction>& cost,
                  const std::shared_ptr<ImpulseConstraints>& constraints);

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
  /// @param[in] impulse_status Impulse status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void linearizeOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                    const double t, const Eigen::VectorXd& q_prev, 
                    const ImpulseSplitSolution& s, 
                    const SplitSolution& s_next);

  template <typename MatrixType, typename VectorType>
  void getStateConstraintFactorization(
      const Eigen::MatrixBase<MatrixType>& T,
      const Eigen::MatrixBase<VectorType>& e) const;

  template <typename MatrixType1, typename MatrixType2>
  void backwardStateConstraintFactorization(
      const Eigen::MatrixBase<MatrixType1>& T_next,
      const Eigen::MatrixBase<MatrixType2>& T) const;

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
  /// @param[in] dtau Length of the discretization of the horizon.
  /// @param[in] d Split direction of this stage.
  /// @param[out] d_next Split direction of the next stage.
  /// 
  void forwardRiccatiRecursionSerial(const RiccatiFactorization& riccati,
                                     RiccatiFactorization& riccati_next);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] riccati Riccati factorization of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// 
  template <typename VectorType>
  void computePrimalDirection(Robot& robot, const RiccatiFactorization& riccati, 
                              const ImpulseSplitSolution& s, 
                              const Eigen::MatrixBase<VectorType>& dx0, 
                              ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeDualDirection(Robot& robot, const SplitDirection& d_next, 
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
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size of the OCP. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] primal_step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(Robot& robot, const double primal_step_size, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  ///
  void computeKKTResidual(Robot& robot, const ImpulseStatus& impulse_status,
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

  ///
  /// @brief Computes the stage cost of this stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] primal_step_size Primal step size of the OCP. Default is 0.
  /// @return Stage cost of this stage.
  /// 
  double stageCost(Robot& robot, const double t, const ImpulseSplitSolution& s, 
                   const double primal_step_size=0);

  ///
  /// @brief Computes and returns the constraint violation of the OCP at this 
  /// stage for line search.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] s Split solution of this stage.
  /// @param[in] q_next Configuration at the next stage.
  /// @param[in] v_next Generaized velocity at the next stage.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ImpulseStatus& impulse_status, 
                             const double t, const ImpulseSplitSolution& s, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<ImpulseCostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<ImpulseConstraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseKKTResidual kkt_residual_;
  ImpulseKKTMatrix kkt_matrix_;
  ImpulseDynamicsForwardEuler impulse_dynamics_;
  ImpulseRiccatiFactorizer riccati_factorizer_;

  ///
  /// @brief Set impulse status from robot model, i.e., set dimension of the 
  /// impulses and equality constraints.
  /// @param[in] impulse_status Impulse status.
  ///
  inline void setImpulseStatusForKKT(const ImpulseStatus& impulse_status) {
    kkt_residual_.setImpulseStatus(impulse_status);
    kkt_matrix_.setImpulseStatus(impulse_status);
  }

};

} // namespace idocp

#include "idocp/impulse/split_impulse_ocp.hxx"

#endif // IDOCP_SPLIT_IMPULSE_OCP_HPP_ 