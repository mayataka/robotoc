#ifndef IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_ 
#define IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class SplitParNMPC
/// @brief Split optimal control problem of ParNMPC algorithm for single stage. 
///
class SplitImpulseParNMPC {
public:
  ///
  /// @brief Construct SplitParNMPC.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  SplitImpulseParNMPC(const Robot& robot, 
                      const std::shared_ptr<CostFunction>& cost,
                      const std::shared_ptr<ImpulseConstraints>& constraints);

  ///
  /// @brief Default constructor. 
  ///
  SplitImpulseParNMPC();

  ///
  /// @brief Destructor. 
  ///
  ~SplitImpulseParNMPC();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitImpulseParNMPC(const SplitImpulseParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitImpulseParNMPC& operator=(const SplitImpulseParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitImpulseParNMPC(SplitImpulseParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitImpulseParNMPC& operator=(SplitImpulseParNMPC&&) noexcept = default;

  ///
  /// @brief Check whether the solution is feasible under inequality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  ///
  bool isFeasible(Robot& robot, const ImpulseSplitSolution& s);

  ///
  /// @brief Initialize the constraints, i.e., set slack and dual variables. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] time_step Time step of this stage.
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
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  /// @param[in] is_state_constraint_valid Specify wheather the pure-state
  /// equality constraint is valid or not.
  ///
  void linearizeOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                    const double t, const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, 
                    const ImpulseSplitSolution& s, 
                    const SplitSolution& s_next, 
                    ImpulseSplitKKTMatrix& kkt_matrix, 
                    ImpulseSplitKKTResidual& kkt_residual,
                    const bool is_state_constraint_valid);

  // ///
  // /// @brief Gets the auxiliary matrix of this stage.
  // /// @param[out] aux_mat The auxiliary matrix of this stage. Size must be 
  // /// 2 * Robot::dimv() x 2 * Robot::dimv().
  // ///
  // void getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const;

  // ///
  // /// @brief Corrects the part of the solution updated coarsely. Call serially 
  // /// after calling coarseUpdate() or coarseUpdateTerminal() and before calling
  // /// backwardCorrectionParallel().
  // /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  // /// @param[in] s_next Split solution of the next stage at the previous 
  // /// iteration.
  // /// @param[in] s_new_next Split solution of the next stage at the current 
  // /// iteration.
  // /// @param[out] s_new Split solution of the current stage at the current 
  // /// iteration.
  // ///
  // void backwardCorrectionSerial(const Robot& robot, const SplitSolution& s_next,
  //                               const SplitSolution& s_next_new,
  //                               ImpulseSplitSolution& s_new);

  // ///
  // /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  // /// after backwardCorrectionSerial() and before calling 
  // /// forwardCorrectionSerial().
  // /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  // /// @param[in] d Split directtion of the current stage.
  // /// @param[out] s_new Split solution of the current stage at the current 
  // /// iteration.
  // ///
  // void backwardCorrectionParallel(const Robot& robot, ImpulseSplitDirection& d,
  //                                 ImpulseSplitSolution& s_new) const;

  // ///
  // /// @brief Corrects the part of the solutio updated coarsely. Call serially 
  // /// after backwardCorrectionParallel() and before calling 
  // /// forwardCorrectionParallel().
  // /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  // /// @param[in] s_prev Split solution of the previous stage at the previous 
  // /// iteration.
  // /// @param[in] s_new_prev Split solution of the previous stage at the current 
  // /// iteration.
  // /// @param[out] s_new Split solution of the current stage at the current 
  // /// iteration.
  // ///
  // void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev,
  //                              const SplitSolution& s_prev_new, 
  //                              ImpulseSplitSolution& s_new);

  // ///
  // /// @brief Corrects the part of the solutio updated coarsely. Call parallel 
  // /// after forwardCorrectionSerial() and before calling 
  // /// computePrimalAndDualDirection().
  // /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  // /// @param[in] d Split directtion of the current stage.
  // /// @param[out] s_new Split solution of the current stage at the current 
  // /// iteration.
  // ///
  // void forwardCorrectionParallel(const Robot& robot, ImpulseSplitDirection& d, 
  //                                ImpulseSplitSolution& s_new) const;

  // ///
  // /// @brief Computes the direction of the primal and dual solution. Call after 
  // /// forwardCorrectionParallel() and before calling 
  // /// computePrimalAndDualDirection().
  // /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  // /// @param[in] s Split solution of this stage.
  // /// @param[in] s_new Corrected split solution of the current stage.
  // /// @param[out] d Split direction of this stage.
  // /// 
  // void computePrimalAndDualDirection(Robot& robot, 
  //                                    const ImpulseSplitSolution& s,
  //                                    const ImpulseSplitSolution& s_new,
  //                                    ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed primal variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] s Split solution of this stage.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedPrimalDirection(Robot& robot,  
                                       const ImpulseSplitSolution& s, 
                                       ImpulseSplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the condensed dual variables of 
  /// this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] kkt_matrix KKT matrix of this stage.
  /// @param[in] kkt_residual KKT residual of this stage.
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in] d Split direction of this stage.
  /// 
  void computeCondensedDualDirection(
      const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
      const ImpulseSplitKKTResidual& kkt_residual, ImpulseSplitDirection& d);

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
  /// @param[in] is_state_constraint_valid Specify wheather the pure-state
  /// equality constraint is valid or not.
  ///
  void updatePrimal(const Robot& robot, const double primal_step_size, 
                    const ImpulseSplitDirection& d, ImpulseSplitSolution& s,
                    const bool is_state_constraint_valid);

  ///
  /// @brief Computes the KKT residual of the OCP at this stage.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status of robot at this stage. 
  /// @param[in] t Current time of this stage. 
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[in] s_next Split solution of the next stage.
  /// @param[out] kkt_matrix KKT matrix of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  /// @param[in] is_state_constraint_valid Specify wheather the pure-state
  /// equality constraint is valid or not.
  ///
  void computeKKTResidual(Robot& robot, const ImpulseStatus& impulse_status, 
                          const double t, const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, 
                          const ImpulseSplitSolution& s, 
                          const SplitSolution& s_next, 
                          ImpulseSplitKKTMatrix& kkt_matrix, 
                          ImpulseSplitKKTResidual& kkt_residual,
                          const bool is_state_constraint_valid);

  ///
  /// @brief Returns the KKT residual of the OCP at this stage. Before calling 
  /// this function, SplitOCP::linearizeOCP or SplitOCP::computeKKTResidual
  /// must be called.
  /// @param[in] kkt_residual KKT residual of this stage.
  /// @param[in] is_state_constraint_valid Specify wheather the pure-state
  /// equality constraint is valid or not.
  /// @return The squared norm of the kKT residual.
  ///
  double squaredNormKKTResidual(const ImpulseSplitKKTResidual& kkt_residual,
                                const bool is_state_constraint_valid) const;

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
  /// @param[in] q_prev Configuration of the previous stage.
  /// @param[in] v_prev Velocity of the previous stage.
  /// @param[in] s Split solution of this stage.
  /// @param[out] kkt_residual KKT residual of this stage.
  /// @param[in] is_state_constraint_valid Specify wheather the pure-state
  /// equality constraint is valid or not.
  /// @return Constraint violation of this stage.
  ///
  double constraintViolation(Robot& robot, const ImpulseStatus& impulse_status, 
                             const double t, const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const ImpulseSplitSolution& s, 
                             ImpulseSplitKKTResidual& kkt_residual,
                             const bool is_state_constraint_valid);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<ImpulseConstraints> constraints_;
  ConstraintsData constraints_data_;
  ImpulseDynamicsBackwardEuler impulse_dynamics_;

};

} // namespace idocp

#include "idocp/impulse/split_impulse_parnmpc.hxx"

#endif // IDOCP_SPLIT_IMPULSE_PARNMPC_HPP_ 