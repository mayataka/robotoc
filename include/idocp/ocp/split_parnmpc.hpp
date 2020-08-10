#ifndef IDParNMPC_SPLIT_PARNMPC_HPP_
#define IDParNMPC_SPLIT_PARNMPC_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class SplitParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  SplitParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
               const std::shared_ptr<ConstraintsInterface>& constraints);

  // Default constructor.
  SplitParNMPC();

  // Destructor.
  ~SplitParNMPC();

  // Use default copy constructor.
  SplitParNMPC(const SplitOCP&) = default;

  // Use default copy assign operator.
  SplitParNMPC& operator=(const SplitOCP&) = default;

  // Use default move constructor.
  SplitParNMPC(SplitOCP&&) noexcept = default;

  // Use default move assign operator.
  SplitParNMPC& operator=(SplitOCP&&) noexcept = default;
 
  // Check whether the solution q, v, a, u are feasible under inequality 
  // constraints.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  bool isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  const Eigen::VectorXd& u);

  // Initialize the constraints, i.e., set slack and dual variables under set 
  //  q, v, a, u.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u);

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
                    const Eigen::VectorXd& v_prev,
                    const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                    const Eigen::VectorXd& mu, const Eigen::VectorXd& a,
                    const Eigen::VectorXd& f, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, const Eigen::VectorXd& u, 
                    const Eigen::VectorXd& lmd_next,
                    const Eigen::VectorXd& gmm_next,
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next,
                    const Eigen::MatrixXd& aux_mat_next_old,
                    Eigen::MatrixXd& aux_mat);

  void backwardCollectionSerial(const Eigen::VectorXd& lmd_next_old, 
                                const Eigen::VectorXd& gmm_next_old, 
                                const Eigen::VectorXd& lmd_next, 
                                const Eigen::VectorXd& gmm_next,
                                Eigen::VectorXd& lmd, 
                                Eigen::VectorXd& gmm);

  void backwardCollectionParallel(const Robot& robot);

  void forwardCollectionSerial(const Robot& robot, const Eigen::VectorXd& q_prev_old,   
                               const Eigen::VectorXd& v_prev_old, 
                               const Eigen::VectorXd& q_prev, 
                               const Eigen::VectorXd& v_prev,
                               Eigen::VectorXd& q, Eigen::VectorXd& v);

  void forwardCollectionParallel(const Robot& robot);

  void computeDirection(const Robot& robot, 
                        const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                        const Eigen::VectorXd& mu, const Eigen::VectorXd& a,
                        const Eigen::VectorXd& f, const Eigen::VectorXd& q, 
                        const Eigen::VectorXd& v, const Eigen::VectorXd& u);

  double maxPrimalStepSize();

  double maxDualStepSize();

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double t, const double dtau, const Eigen::VectorXd& q, 
      const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
      const Eigen::VectorXd& f, const Eigen::VectorXd& u);

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& f, 
      const Eigen::VectorXd& u, const Eigen::VectorXd& q_next, 
      const Eigen::VectorXd& v_next, const Eigen::VectorXd& dq, 
      const Eigen::VectorXd& dv, const Eigen::VectorXd& dq_next, 
      const Eigen::VectorXd& dv_next);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a, 
                    Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                    Eigen::VectorXd& f, Eigen::VectorXd& mu);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& q_prev, 
                             const Eigen::VectorXd& v_prev,
                             const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& mu, const Eigen::VectorXd& a,
                             const Eigen::VectorXd& f, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, const Eigen::VectorXd& u, 
                             const Eigen::VectorXd& lmd_next,
                             const Eigen::VectorXd& gmm_next,
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

private:
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<ConstraintsInterface> constraints_;
  CostFunctionData cost_data_;
  pdipm::JointSpaceConstraints joint_constraints_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  KKTComposition kkt_composition_;
  bool has_floating_base_;
  int dimq_, dimv_, dim_passive_, max_dimf_, max_dimc_;
  Eigen::VectorXd lu_, lu_condensed_, u_res_, du_, x_res_, dx_, 
                  dkkt_, dlmd_, dgmm_, dmu_, da_, df_, dq_, dv_, du_, dbeta_;
  Eigen::MatrixXd luu_, kkt_matrix_inverse_, du_dq_, du_dv_, du_da_, du_df_,
                  dsubtract_dqminus_, dsubtract_dqplus_;
  Eigen::VectorXd  lmd_tmp_, gmm_tmp_, q_tmp_, v_tmp_, a_tmp_, f_tmp_, u_tmp_, u_res_tmp_;
};

} // namespace idocp


#endif // IDParNMPC_SPLIT_PARNMPC_HPP_