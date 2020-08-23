#ifndef IDOCP_SPLIT_PARNMPC_HPP_ 
#define IDOCP_SPLIT_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/robot_dynamics.hpp"


namespace idocp {

class SplitParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The shared pointer to the cost function.
  //    constraints: The shared pointer to the inequality constraints.
  SplitParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
               const std::shared_ptr<Constraints>& constraints);

  // Default constructor.
  SplitParNMPC();

  // Destructor.
  ~SplitParNMPC();

  // Use default copy constructor.
  SplitParNMPC(const SplitParNMPC&) = default;

  // Use default copy assign operator.
  SplitParNMPC& operator=(const SplitParNMPC&) = default;

  // Use default move constructor.
  SplitParNMPC(SplitParNMPC&&) noexcept = default;

  // Use default move assign operator.
  SplitParNMPC& operator=(SplitParNMPC&&) noexcept = default;
 
  // Check whether the solution s is feasible under inequality constraints.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   s: split solution of this stage.
  bool isFeasible(const Robot& robot, const SplitSolution& s);

  // Initialize the constraints, i.e., set slack and dual variables.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   time_step: The current discreization time step.
  //   dtau: The length of the discretization.
  //   s: split solution of this stage.
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  // Updates the solution approximately.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   q_prev: Configuration at the previous stage. Size must be robot.dimq().
  //   v_prev: Generalized velocity at the previous stage. Size must be 
  //      robot.dimv().
  //   s: Split solution of this stage.
  //   s_next: Split solution of the next stage.
  //   aux_mat_next_old: Guess of the auxiliary matrix of the next stage. The 
  //      size must be 2*robot.dimv() x 2*robot.dimv().
  //   d: Split direction of this stage. 
  //   s_new_coarse: Coarse updated split solution of this stage.
  void coarseUpdate(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                    const SplitSolution& s_next, 
                    const Eigen::MatrixXd& aux_mat_next_old, SplitDirection& d, 
                    SplitSolution& s_new_coarse);

  // Updates the solution approximately at the terminal stage.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   q_prev: Configuration at the previous stage. Size must be robot.dimq().
  //   v_prev: Generalized velocity at the previous stage. Size must be 
  //      robot.dimv().
  //   s: Split solution of this stage.
  //   d: Split direction of this stage. 
  //   s_new_coarse: Coarse updated split solution of this stage.
  void coarseUpdateTerminal(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& q_prev, 
                            const Eigen::VectorXd& v_prev, 
                            const SplitSolution& s, SplitDirection& d, 
                            SplitSolution& s_new_coarse);

  // Getter of the auxiliary matrix at this stage.
  // Argments: 
  //   aux_mat: The auxiliary matrix of this stage. The size must be 
  //      2*dimv() x 2*dimv().
  void getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const;

  // Getter of the terminal cost. Useful to initialize the auxiliary matrix.
  // Argments: 
  //   robot: The robot model. 
  //   t: Time of the current stage.
  //   s: Split solution of this stage.
  //   phixx: The hessian of the terminal cost. The size must be 
  //      2*robot.dimv() x 2*robot.dimv().
  void getTerminalCostHessian(Robot& robot, const double t, 
                              const SplitSolution& s, Eigen::MatrixXd& phixx);

  // Correct the part of the solutio updated coarsely. Call serially after 
  // calling coarseUpdate() or coarseUpdateTerminal() and before calling
  // backwardCorrectionParallel().
  // Argments: 
  //   robot: The robot model. 
  //   s_next: Split solution of the next stage at the previous iteration.
  //   s_new_next: Split solution of the next stage at the current iteration.
  //   s_new: Split solution of the current stage at the current iteration.
  void backwardCorrectionSerial(const Robot& robot, const SplitSolution& s_next,
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  // Correct the part of the solutio updated coarsely. Call parallel after 
  // backwardCorrectionSerial() and before calling forwardCorrectionSerial().
  // Argments: 
  //   robot: The robot model. 
  //   d: Split directtion of the current stage.
  //   s_new: Split solution of the current stage at the current iteration.
  void backwardCorrectionParallel(const Robot& robot, SplitDirection& d,
                                  SplitSolution& s_new);

  // Correct the part of the solutio updated coarsely. Call serially after 
  // backwardCorrectionParallel() and before calling forwardCorrectionParallel().
  // Argments: 
  //   robot: The robot model. 
  //   s_prev: Split solution of the previous stage at the previous iteration.
  //   s_new_prev: Split solution of the previous stage at the current iteration.
  //   s_new: Split solution of the current stage at the current iteration.
  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev,
                               const SplitSolution& s_new_prev, 
                               SplitSolution& s_new);

  // Correct the part of the solutio updated coarsely. Call parallel after 
  // forwardCorrectionSerial() and before calling computePrimalAndDualDirection().
  // Argments: 
  //   robot: The robot model. 
  //   d: Split directtion of the current stage.
  //   s_new: Split solution of the current stage at the current iteration.
  void forwardCorrectionParallel(const Robot& robot, SplitDirection& d, 
                                 SplitSolution& s_new);

  // Computes the direction of the primal and dual solution. Call after 
  // forwardCorrectionParallel() and before calling 
  // computePrimalAndDualDirection().
  // Argments: 
  //   robot: The robot model. 
  //   dtau: Discretization length of the OCP. Must be positive.
  //   s: Split solution of this stage at the previous iteration.
  //   s_new: Corrected split solution of the current stage.
  //   d: Computed split directtion of the current stage.
  void computePrimalAndDualDirection(const Robot& robot, const double dtau,
                                     const SplitSolution& s,
                                     const SplitSolution& s_new,
                                     SplitDirection& d);

  // Returns the max step size of the primal solution.
  double maxPrimalStepSize();

  // Returns the max step size of the dual solution.
  double maxDualStepSize();

  // Returns the stage cost and the l1 norm of the constraints violation at 
  //    this stage.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   s: Split solution of this stage.
  std::pair<double, double> costAndViolation(Robot& robot, const double t, 
                                             const double dtau, 
                                             const SplitSolution& s);

  // Returns the stage cost and the l1 norm of the constraints violation at 
  //    this stage under the step size.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   step_size: Step size candidate of the solution update.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   q_prev: Configuration at the previous stage. Size must be robot.dimq().
  //   v_prev: Generalized velocity at the previous stage. Size must be 
  //      robot.dimv().
  //   s: Split solution of this stage.
  //   d: Split direction of this stage. 
  //   s_tmp: Tempolary split solution.
  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const Eigen::VectorXd& q_prev, 
                                             const Eigen::VectorXd& v_prev, 
                                             const SplitSolution& s, 
                                             const SplitDirection& d, 
                                             SplitSolution& s_tmp);

  // Returns the stage cost and the l1 norm of the constraints violation at 
  //    this stage under the step size.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   step_size: Step size candidate of the solution update.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   s_prev: Split solution of the previous stage.
  //   d_prev: Split direction of the previous stage.
  //   s: Split solution of this stage.
  //   d: Split direction of this stage. 
  //   s_tmp: Tempolary split solution.
  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const SplitSolution& s_prev, 
                                             const SplitDirection& d_prev, 
                                             const SplitSolution& s, 
                                             const SplitDirection& d, 
                                             SplitSolution& s_tmp);

  // Returns the stage cost and the l1 norm of the constraints violation at 
  //    this terminal stage.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   s: Split solution of this stage.
  std::pair<double, double> costAndViolationTerminal(Robot& robot, 
                                                     const double t, 
                                                     const double dtau, 
                                                     const SplitSolution& s);

  // Returns the stage cost and the l1 norm of the constraints violation at 
  //    this stage under the step size.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   step_size: Step size candidate of the solution update.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   s_prev: Split solution of the previous stage.
  //   d_prev: Split direction of the previous stage.
  //   s: Split solution of this stage.
  //   d: Split direction of this stage. 
  //   s_tmp: Tempolary split solution.
  std::pair<double, double> costAndViolationTerminal(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const SplitSolution& s_prev, const SplitDirection& d_prev,
      const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp);

  // Updates the primal solution.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   step_size: Step size of the solution update.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   d: Split direction of this stage. 
  //   s: Split solution of this stage.
  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const SplitDirection& d, SplitSolution& s);

  // Updates the dual solution.
  // Argments: 
  //   step_size: Step size of the solution update.
  void updateDual(const double step_size);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  // Returns the squared KKT error norm of the current stage of the split OCP.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   q_prev: Configuration at the previous stage. Size must be robot.dimq().
  //   v_prev: Generalized velocity at the previous stage. Size must be 
  //      robot.dimv().
  //   s: Split solution of this stage.
  //   s_next: Split solution of the next stage.
  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const SplitSolution& s,
                             const SplitSolution& s_next);

  // Returns the squared KKT error norm of the current terminal stage of the 
  // split OCP.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      must be included in this model.
  //   t: Time of the current stage.
  //   dtau: Discretization length of the OCP. Must be positive.
  //   q_prev: Configuration at the previous stage. Size must be robot.dimq().
  //   v_prev: Generalized velocity at the previous stage. Size must be 
  //      robot.dimv().
  //   s: split Solution of this stage.
  double squaredKKTErrorNormTerminal(Robot& robot, const double t, 
                                     const double dtau, 
                                     const Eigen::VectorXd& q_prev, 
                                     const Eigen::VectorXd& v_prev, 
                                     const SplitSolution& s);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  StateEquation state_equation_;
  RobotDynamics robot_dynamics_;
  int dimx_, dimKKT_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  Eigen::VectorXd x_res_, dx_;

};

} // namespace idocp


#endif // IDOCP_SPLIT_PARNMPC_HPP_