#ifndef IDParNMPC_SPLIT_PARNMPC_HPP_
#define IDParNMPC_SPLIT_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_direction.hpp"
#include "idocp/ocp/kkt_composition.hpp"
#include "idocp/ocp/split_solution.hpp"


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
                    const Eigen::VectorXd& lmd_next,
                    const Eigen::VectorXd& gmm_next,
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next,
                    const Eigen::MatrixXd& aux_mat_next_old,
                    Eigen::MatrixXd& aux_mat, SplitSolution& s_new_coarse, 
                    const bool is_terminal_ocp=false);

  void computeTerminalCostDerivatives(const Robot& robot, const double t, 
                                      const SplitSolution& s,
                                      Eigen::VectorXd& phiq, 
                                      Eigen::VectorXd& phiv);


  void backwardCollectionSerial(const SplitSolution& s_old_next,
                                const SplitSolution& s_new_next,
                                SplitSolution& s_new);

  void backwardCollectionParallel(const Robot& robot, SplitSolution& s_new);

  void forwardCollectionSerial(const Robot& robot, 
                               const SplitSolution& s_old_prev,
                               const SplitSolution& s_new_prev, 
                               SplitSolution& s_new);

  void forwardCollectionParallel(const Robot& robot, SplitSolution& s_new);

  void computePrimalAndDualDirection(const Robot& robot, const double dtau,
                                     const SplitSolution& s,
                                     const SplitSolution& s_new);

  double maxPrimalStepSize();

  double maxDualStepSize();

  // std::pair<double, double> stageCostAndConstraintsViolation(
  //     Robot& robot, const double t, const double dtau, const SplitSolution& s);

  // std::pair<double, double> stageCostAndConstraintsViolation(
  //     Robot& robot, const double step_size, const double t, const double dtau, 
  //     const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
  //     const Eigen::VectorXd& dq_prev, const Eigen::VectorXd& dv_prev, 
  //     const SplitSolution& s, SplitSolution& s_tmp);

  // double terminalCost(Robot& robot, const double t, const SplitSolution& s);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    SplitSolution& s);

  void updateDual(const double step_size);

  void getStateDirection(Eigen::VectorXd& dq, Eigen::VectorXd& dv);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev, 
                             const SplitSolution& s,
                             const Eigen::VectorXd& lmd_next,
                             const Eigen::VectorXd& gmm_next,
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);
  

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  pdipm::JointSpaceConstraints joint_constraints_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  KKTDirection kkt_direction_;
  KKTComposition kkt_composition_;
  Eigen::VectorXd lu_, lu_condensed_, u_res_, du_, dbeta_, x_res_, dx_,
                  u_tmp_, u_res_tmp_;
  Eigen::MatrixXd luu_, kkt_matrix_inverse_, du_dq_, du_dv_, du_da_, du_df_,
                  dsubtract_dqminus_, dsubtract_dqplus_;
};

} // namespace idocp


#endif // IDParNMPC_SPLIT_PARNMPC_HPP_