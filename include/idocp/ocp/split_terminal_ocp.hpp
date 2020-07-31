#ifndef IDOCP_SPLIT_TERMINAL_OCP_HPP_
#define IDOCP_SPLIT_TERMINAL_OCP_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "constraints/joint_space_constraints/joint_space_constraints.hpp"


namespace idocp {

class SplitTerminalOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    cost: The pointer to the constraints.
  SplitTerminalOCP(const Robot& robot, const CostFunctionInterface* cost,
                   const ConstraintsInterface* constraints);

  // Destructor.
  ~SplitTerminalOCP();
 
  // Use default copy constructor.
  SplitTerminalOCP(const SplitTerminalOCP& other) = default;

  // Use default copy operator.
  SplitTerminalOCP& operator=(const SplitTerminalOCP& other) = default;

  // Check whether the solution q, v, a, u, fext are feasible under inequality 
  // constraints.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  bool isFeasible(Robot& robot, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  // Initialize the constraints, i.e., set slack and dual variables under set 
  //  q, v, a, u.
  // Argments: 
  //   robot: The robot model that has been already initialized.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  void initConstraints(Robot& robot, const double dtau,
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v);

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

  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv);
 
  double maxPrimalStepSize();

  double maxDualStepSize();

  double costDerivativesDotDirection(Robot& robot, const double t, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const Eigen::VectorXd& dq,
                                     const Eigen::VectorXd& dv);

  double terminalCost(Robot& robot, const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& dq, const Eigen::VectorXd& dv);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v);

private:
  CostFunctionInterface *cost_;
  ConstraintsInterface *constraints_;
  pdipm::JointSpaceConstraints joint_constraints_;
  int dimq_, dimv_;
  Eigen::VectorXd lq_, lv_, q_res_, v_res_;
  // The following variables are only needed for line search
  Eigen::VectorXd q_tmp_, v_tmp_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_TERMINAL_OCP_HPP_