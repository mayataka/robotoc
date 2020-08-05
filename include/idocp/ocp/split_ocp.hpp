#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

#include <utility>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class SplitOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  SplitOCP(const Robot& robot, 
           std::unique_ptr<CostFunctionInterface> cost,
           std::unique_ptr<ConstraintsInterface> constraints);

  // Default constructor.
  SplitOCP();

  // Destructor.
  ~SplitOCP();

  // Progibits copy constructor due to unique_ptr.
  SplitOCP(const SplitOCP&) = delete;

  // Progibits copy operator due to unique_ptr.
  SplitOCP& operator=(const SplitOCP&) = delete;

  // Use default move constructor.
  SplitOCP(SplitOCP&&) noexcept = default;

  // Use default move assign operator.
  SplitOCP& operator=(SplitOCP&&) noexcept = default;
 
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

  // Linearize the OCP for Newton's method around the current solution.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   dtau: Discretization length of the OCP.
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
  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& lmd, const Eigen::VectorXd& gmm, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                    const Eigen::VectorXd& f, const Eigen::VectorXd& mu,
                    const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& gmm_next, 
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next);

  void backwardRiccatiRecursion(const double dtau, 
                                const Eigen::MatrixXd& Pqq_next, 
                                const Eigen::MatrixXd& Pqv_next, 
                                const Eigen::MatrixXd& Pvq_next, 
                                const Eigen::MatrixXd& Pvv_next, 
                                const Eigen::VectorXd& sq_next, 
                                const Eigen::VectorXd& sv_next, 
                                Eigen::MatrixXd& Pqq, Eigen::MatrixXd& Pqv, 
                                Eigen::MatrixXd& Pvq, Eigen::MatrixXd& Pvv, 
                                Eigen::VectorXd& sq, Eigen::VectorXd& sv);

  void forwardRiccatiRecursion(const double dtau, const Eigen::VectorXd& dq,   
                               const Eigen::VectorXd& dv, 
                               Eigen::VectorXd& dq_next, 
                               Eigen::VectorXd& dv_next);

  void computeCondensedDirection(const double dtau, const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv);
 
  double maxPrimalStepSize();

  double maxDualStepSize();

  std::pair<double, double> costAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
      const Eigen::VectorXd& f, const Eigen::VectorXd& q_next, 
      const Eigen::VectorXd& v_next, const Eigen::VectorXd& dq, 
      const Eigen::VectorXd& dv, const Eigen::VectorXd& dq_next, 
      const Eigen::VectorXd& dv_next);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a, 
                    Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                    Eigen::VectorXd& f, Eigen::VectorXd& mu, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                             const Eigen::VectorXd& beta, 
                             const Eigen::VectorXd& f, const Eigen::VectorXd& mu, 
                             const Eigen::VectorXd& lmd_next, 
                             const Eigen::VectorXd& gmm_next, 
                             const Eigen::VectorXd& q_next, 
                             const Eigen::VectorXd& v_next);

private:
  std::unique_ptr<CostFunctionInterface> cost_;
  std::unique_ptr<ConstraintsInterface> constraints_;
  pdipm::JointSpaceConstraints joint_constraints_;
  RiccatiMatrixFactorizer riccati_matrix_factorizer_;
  RiccatiMatrixInverter riccati_matrix_inverter_;
  bool has_floating_base_;
  int dimq_, dimv_, dim_passive_, max_dimf_, max_dimc_, dimf_, dimc_;
  Eigen::VectorXd lq_, lv_, la_, lf_, lu_, lu_condensed_, ka_, kf_,  kmu_, 
                  da_, df_, dmu_, q_res_, v_res_, u_res_, du_, C_res_;
  Eigen::MatrixXd luu_, du_dq_, du_dv_, du_da_, du_df_, Qqq_, Qqv_, Qqa_, Qqf_, 
                  Qvq_, Qvv_, Qva_, Qvf_, Qaa_, Qaf_, Qff_, Cq_, Cv_, Ca_, Cf_, 
                  Kaq_, Kav_, Kfq_, Kfv_, Kmuq_, Kmuv_;
  // The following variables are only needed for line search
  Eigen::VectorXd q_tmp_, v_tmp_, a_tmp_, f_tmp_, u_tmp_, u_res_tmp_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_