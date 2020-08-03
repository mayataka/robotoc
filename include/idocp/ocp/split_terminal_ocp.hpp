#ifndef IDOCP_SPLIT_TERMINAL_OCP_HPP_
#define IDOCP_SPLIT_TERMINAL_OCP_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_interface.hpp"
#include "idocp/constraints/constraints_interface.hpp"
#include "idocp/constraints/joint_space_constraints/joint_space_constraints.hpp"


namespace idocp {

class SplitTerminalOCP {
public:
  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  SplitTerminalOCP(const Robot& robot, const CostFunctionInterface* cost,
                   const ConstraintsInterface* constraints);

  // Default constructor. 
  SplitTerminalOCP();

  // Destructor.
  ~SplitTerminalOCP();
 
  // Use default copy constructor.
  SplitTerminalOCP(const SplitTerminalOCP& other) = default;

  // Use default copy operator.
  SplitTerminalOCP& operator=(const SplitTerminalOCP& other) = default;

  // Sets the cost function.
  // Argments:
  //    cost: The pointer to the cost function.
  void setCostFunction(const CostFunctionInterface* cost);

  // Sets the constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    constraints: The pointer to the constraints.
  void setConstraints(const ConstraintsInterface* constraints);

  // Check whether the solution q, v are feasible under inequality 
  // constraints.
  // Argments: 
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  bool isFeasible(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  // Initialize the constraints, i.e., set slack and dual variables under set 
  //  q, v.
  // Argments: 
  //   time_step: The time step of the split OCP.
  //   dtau: Discretization interval of the horizon.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  void initConstraints(const int time_step, const double dtau,
                       const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  // Linearize the OCP for Newton's method around the current solution at the 
  // last time step of the horizon.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   dtau: Discretization interval of the horizon.
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

  // Computes the direction of the condensed veriables.
  // Argments: 
  //   dtau: Discretization interval of the horizon.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv);

  // Returns the maximum step size of the primal variables of the inequality 
  // constraints.
  double maxPrimalStepSize();

  // Returns the maximum step size of the dual variables of the inequality 
  // constraints.
  double maxDualStepSize();

  // Returns the terminal cost.
  // Argments: 
  //   t: Time of the current time step.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  double terminalCost(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  // Returns the terminal cost with step_size.
  // Argments: 
  //   robot: The robot model. 
  //   step_size: The step size.
  //   t: Time of the current time step.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                      const Eigen::VectorXd& dq, const Eigen::VectorXd& dv);

  // Updates the dual variables of the inequality constraints.
  // Argments: 
  //   step_size: The step size.
  void updateDual(const double step_size);

  // Updates the primal variables.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   step_size: The step size.
  //   Pqq: The Riccati factorization matrix. The size must be dimv x dimv.
  //   Pqv: The Riccati factorization matrix. The size must be dimv x dimv.
  //   Pvq: The Riccati factorization matrix. The size must be dimv x dimv.
  //   Pvv: The Riccati factorization matrix. The size must be dimv x dimv.
  //   sq: The Riccati factorization vector. The size must be dimv.
  //   sv: The Riccati factorization vector. The size must be dimv.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  //   lmd: The Lagrange multiplier with respect to the transition of the 
  //      configuration. Size must be dimv.
  //   gmm: The Lagrange multiplier with respect to the transition of the 
  //      generalized velocity. Size must be dimv.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  void updatePrimal(Robot& robot, const double step_size, 
                    const Eigen::MatrixXd& Pqq, const Eigen::MatrixXd& Pqv, 
                    const Eigen::MatrixXd& Pvq, const Eigen::MatrixXd& Pvv, 
                    const Eigen::VectorXd& sq, const Eigen::VectorXd& sv, 
                    const Eigen::VectorXd& dq, const Eigen::VectorXd& dv, 
                    Eigen::VectorXd& lmd, Eigen::VectorXd& gmm, 
                    Eigen::VectorXd& q, Eigen::VectorXd& v) const;

  // Returns the squared KKT error norm.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   lmd: The Lagrange multiplier with respect to the transition of the 
  //      configuration. Size must be dimv.
  //   gmm: The Lagrange multiplier with respect to the transition of the 
  //      generalized velocity. Size must be dimv.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  double squaredKKTErrorNorm(Robot& robot, const double t, 
                             const Eigen::VectorXd& lmd, 
                             const Eigen::VectorXd& gmm, 
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v);

private:
  CostFunctionInterface* cost_;
  ConstraintsInterface* constraints_;
  pdipm::JointSpaceConstraints joint_constraints_;
  int dimq_, dimv_;
  Eigen::VectorXd lq_, lv_, q_res_, v_res_;
  // The following variables are only needed for line search
  Eigen::VectorXd q_tmp_, v_tmp_;
};

} // namespace idocp


#endif // IDOCP_SPLIT_TERMINAL_OCP_HPP_