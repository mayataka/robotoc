#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "ocp/riccati_matrix_factorizer.hpp"
#include "ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class SplitOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    cost: The pointer to the constraints.
  SplitOCP(const Robot& robot, const CostFunctionInterface* cost,
           const ConstraintsInterface* constraints);

  // Destructor.
  ~SplitOCP();
 
  // Use default copy constructor.
  SplitOCP(const SplitOCP& other) = default;

  // Use default copy operator.
  SplitOCP& operator=(const SplitOCP& other) = default;

  // Check whether the solution q, v, a, u, fext are feasible under inequality 
  // constraints.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  bool isFeasible(Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  const Eigen::VectorXd& u);

  // Initialize the constraints, i.e., set slack and dual variables under set 
  //  q, v, a, u.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   u: Generalized torque. Size must be dimv.
  void initConstraints(Robot& robot, const int time_step, const double dtau, 
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                       const Eigen::VectorXd& a, const Eigen::VectorXd& u);

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
                    const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& gmm_next, 
                    const Eigen::VectorXd& q_next,
                    const Eigen::VectorXd& v_next);

  // Linearize the OCP for Newton's method around the current solution at the 
  // last time step of the horizon.
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
  //   Qqq: The second partial derivative of the terminal cost with qq. Size 
  //      must be dimv x dimv.
  //   Qqv: The second partial derivative of the terminal cost with qv. Size
  //      must be dimv x dimv.
  //   Qvq: The second partial derivative of the terminal cost with vq. Size
  //      must be dimv x dimv.
  //   Qvv: The second partial derivative of the terminal cost with vv. Size
  //      must be dimv x dimv.
  //   Qq: The partial derivative of the terminal cost with q. Size must be dimv.
  //   Qv: The partial derivative of the terminal cost with v. Size must be dimv. 
  void linearizeOCP(Robot& robot, const double t, const Eigen::VectorXd& lmd, 
                    const Eigen::VectorXd& gmm, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::MatrixXd& Qqq, 
                    Eigen::MatrixXd& Qqv, Eigen::MatrixXd& Qvq, 
                    Eigen::MatrixXd& Qvv, Eigen::VectorXd& Qq, 
                    Eigen::VectorXd& Qv);

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

  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv);
 
  double maxPrimalStepSize();

  double maxDualStepSize();

  double costDerivativeDotDirection(Robot& robot, const double t, 
                                    const double dtau, const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v, 
                                    const Eigen::VectorXd& a, 
                                    const Eigen::VectorXd& u,
                                    const Eigen::VectorXd& dq,
                                    const Eigen::VectorXd& dv);

  std::pair<double, double> stageCostAndConstraintsViolation(
      Robot& robot, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  std::pair<double, double> stageCostAndConstraintsViolation(
      Robot& robot, const double step_size, const double t, const double dtau, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
      const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
      const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
      const Eigen::VectorXd& dq_next, const Eigen::VectorXd& dv_next);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, Eigen::VectorXd& a, 
                    Eigen::VectorXd& u, Eigen::VectorXd& beta, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                             const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                             const Eigen::VectorXd& beta, 
                             const Eigen::VectorXd& lmd_next, 
                             const Eigen::VectorXd& gmm_next, 
                             const Eigen::VectorXd& q_next,
                             const Eigen::VectorXd& v_next);

private:
  CostFunctionInterface *cost_;
  ConstraintsInterface *constraints_;
  pdipm::JointSpaceConstraints joint_constraints_;
  RiccatiMatrixFactorizer riccati_matrix_factorizer_;
  RiccatiMatrixInverter riccati_matrix_inverter_;
  int dimq_, dimv_, dimf_, dim_passive_;
  Eigen::VectorXd f_, mu_, lq_, lv_, la_, lf_, lu_, lu_condensed_, 
                  ka_, kf_, kmu_, da_, df_, dmu_, q_res_, v_res_, a_res_, 
                  f_res_, u_res_, du_, C_res_;
  Eigen::MatrixXd luu_, du_dq_, du_dv_, du_da_, du_df_, Qqq_, Qqv_, Qqa_, Qqf_, 
                  Qvq_, Qvv_, Qva_, Qvf_, Qaa_, Qaf_, Qff_, Cq_, Cv_, Ca_, Cf_, 
                  Kaq_, Kav_, Kfq_, Kfv_, Kmuq_, Kmuv_;
  // The following variables are only needed for line search
  Eigen::VectorXd q_tmp_, v_tmp_, a_tmp_, f_tmp_, u_tmp_, u_res_tmp_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_