#ifndef IDOCP_IMPULSE_SPLIT_IMPULSE_HPP_ 
#define IDOCP_IMPULSE_SPLIT_IMPULSE_HPP_

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
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class ImpulseSplitParNMPC
/// @brief Split optimal control problem of ParNMPC algorithm for single stage. 
///
class ImpulseSplitParNMPC {
public:
  ///
  /// @brief Construct ImpulseSplitParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  ImpulseSplitParNMPC(const Robot& robot, 
                      const std::shared_ptr<CostFunction>& cost,
                      const std::shared_ptr<Constraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitParNMPC(const ImpulseSplitParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ImpulseSplitParNMPC& operator=(const ImpulseSplitParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitParNMPC(ImpulseSplitParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitParNMPC& operator=(ImpulseSplitParNMPC&&) noexcept = default;

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
  /// @param[in] q_prev Configuration at the previous stage.
  /// @param[in] v_prev Generalized velocity at the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  ///
  void linearizeOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                    const double t, const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, 
                    const ImpulseSplitSolution& s, 
                    const SplitSolution& s_next, 
                    ImpulseSplitKKTMatrix& kkt_matrix, 
                    ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] kkt_matrix KKT matrix of this stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot,  
                                       const ImpulseSplitKKTMatrix& kkt_matrix,
                                       const ImpulseSplitSolution& s, 
                                       ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in, out] d Split direction of this stage.
  /// 
  void computeCondensedDualDirection(const Robot& robot, 
                                     const ImpulseSplitKKTMatrix& kkt_matrix, 
                                     ImpulseSplitKKTResidual& kkt_residual,
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
  /// @brief Updates primal variables of this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] primal_step_size Primal step size of the OCP. 
  /// @param[in] d Split direction of this stage.
  /// @param[in, out] s Split solution of this stage.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s,
                    const bool is_position_constraint_valid);

  ///
  /// @brief Updates dual variables of the inequality constraints.
  /// @param[in] dual_step_size Dula step size of the OCP. 
  ///
  void updateDual(const double dual_step_size);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration at the previous stage.
  /// @param[in] v_prev Generalized velocity at the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  ///
  void computeKKTResidual(Robot& robot, const ImpulseStatus& impulse_status, 
                          const double t, const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, 
                          const ImpulseSplitSolution& s, 
                          const SplitSolution& s_next, 
                          ImpulseSplitKKTMatrix& kkt_matrix, 
                          ImpulseSplitKKTResidual& kkt_residual);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @param[in] kkt_residual KKT residual of this stage.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

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
  /// @param[in] q_prev Configuration at the previous stage.
  /// @param[in] v_prev Generalized velocity at the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ImpulseStatus& impulse_status, 
                             const double t, const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const ImpulseSplitSolution& s, 
                             ImpulseSplitKKTResidual& kkt_residual);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseDynamicsBackwardEuler impulse_dynamics_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_parnmpc.hxx"

#endif // IDOCP_IMPULSE_SPLIT_IMPULSE_HPP_ 