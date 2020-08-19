#ifndef IDParNMPC_SPLIT_PARNMPC_HPP_
#define IDParNMPC_SPLIT_PARNMPC_HPP_

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
#include "idocp/ocp/inverse_dynamics.hpp"
#include "idocp/ocp/equality_constraints.hpp"


namespace idocp {

class SplitParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  //
  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  //
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
 
  // Check whether the solution q, v, a, u are feasible under inequality 
  // constraints.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  bool isFeasible(const Robot& robot, const SplitSolution& s);

  // Initialize the constraints, i.e., set slack and dual variables under set 
  //  q, v, a, u.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  // Linearize the ParNMPC for Newton's method around the current solution.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   dtau: Discretization length of the ParNMPC.
  //   lmd: The Lagrange multiplier with respect to the transition of the 
  //      configuration. Size must be dimv.
  //   gmm: The Lagrange multiplier with respect to the transition of the 
  //      generalized velocity. Size must be dimv.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  //   lmd_next: The Lagrange multiplier with respect to the transition of the 
  //      configuration at the next time step. Size must be dimv.
  //   gmm_next: The Lagrange multiplier with respect to the transition of the 
  //      generalized velocity at the next time step. Size must be dimv.
  //   q_next: Configuration at the next time step. Size must be dimq.
  //   v_next: Generalized velocity at the next time step. Size must be dimv.
  void coarseUpdate(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& q_prev, 
                    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
                    const SplitSolution& s_next, 
                    const Eigen::MatrixXd& aux_mat_next_old, SplitDirection& d, 
                    SplitSolution& s_new_coarse);

  void coarseUpdateTerminal(Robot& robot, const double t, const double dtau, 
                            const Eigen::VectorXd& q_prev, 
                            const Eigen::VectorXd& v_prev, 
                            const SplitSolution& s, SplitDirection& d, 
                            SplitSolution& s_new_coarse);

  void getAuxiliaryMatrix(Eigen::MatrixXd& auxiliary_matrix) const;

  void backwardCorrectionSerial(const Robot& robot, const SplitSolution& s_next,
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  void backwardCorrectionParallel(const Robot& robot, SplitDirection& d,
                                  SplitSolution& s_new);

  void forwardCorrectionSerial(const Robot& robot, const SplitSolution& s_prev,
                               const SplitSolution& s_new_prev, 
                               SplitSolution& s_new);

  void forwardCorrectionParallel(const Robot& robot, SplitDirection& d, 
                                 SplitSolution& s_new);

  void computePrimalAndDualDirection(const Robot& robot, const double dtau,
                                     const SplitSolution& s,
                                     const SplitSolution& s_new,
                                     SplitDirection& d);

  double maxPrimalStepSize();

  double maxDualStepSize();

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double t, const double dtau, const SplitSolution& s);

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
      const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp);

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const SplitSolution& s_prev, const SplitDirection& d_prev,
      const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp);

  std::pair<double, double> costAndConstraintsViolationTerminal(
      Robot& robot, const double t, const double dtau, const SplitSolution& s);

  std::pair<double, double> costAndConstraintsViolationTerminal(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const SplitSolution& s_prev, const SplitDirection& d_prev,
      const SplitSolution& s, const SplitDirection& d, SplitSolution& s_tmp);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const SplitDirection& d, SplitSolution& s);

  void updateDual(const double step_size);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const SplitSolution& s,
                             const SplitSolution& s_next);

  double squaredKKTErrorNormTerminal(Robot& robot, const double t, 
                                     const double dtau, 
                                     const Eigen::VectorXd& q_prev, 
                                     const Eigen::VectorXd& v_prev, 
                                     const SplitSolution& s);

  void setRegularization(const double regularization);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  StateEquation state_equation_;
  InverseDynamics inverse_dynamics_;
  bool use_regularization_;
  double regularization_;
  int dimx_, dimKKT_;
  Eigen::MatrixXd kkt_matrix_inverse_;
  Eigen::VectorXd x_res_, dx_;

  void computeKKTResidual(Robot& robot, const double t, const double dtau, 
                          const Eigen::VectorXd& q_prev, 
                          const Eigen::VectorXd& v_prev, 
                          const SplitSolution& s,
                          const SplitSolution& s_next);

  void computeKKTResidualTerminal(Robot& robot, const double t, 
                                  const double dtau, 
                                  const Eigen::VectorXd& q_prev, 
                                  const Eigen::VectorXd& v_prev, 
                                  const SplitSolution& s);

};

} // namespace idocp


#endif // IDParNMPC_SPLIT_PARNMPC_HPP_