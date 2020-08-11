#ifndef IDOCP_TERMINAL_PARNMPC_HPP_
#define IDOCP_TERMINAL_PARNMPC_HPP_

#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class TerminalParNMPC {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  TerminalParNMPC(const Robot& robot, 
                  const std::shared_ptr<CostFunctionInterface>& cost);

  // Default constructor.
  TerminalParNMPC();
  
  // Destructor.
  ~TerminalParNMPC();

  // Use default copy constructor.
  TerminalParNMPC(const TerminalParNMPC&) = default;

  // Use default copy assign operator.
  TerminalParNMPC& operator=(const TerminalParNMPC&) = default;

  // Use default move constructor.
  TerminalParNMPC(TerminalParNMPC&&) noexcept = default;

  // Use default move assign operator.
  TerminalParNMPC& operator=(TerminalParNMPC&&) noexcept = default;

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
  void linearizeOCP(Robot& robot, const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::VectorXd& phiq, 
                    Eigen::VectorXd& phiv);

  // Returns the terminal cost.
  // Argments: 
  //   t: Time of the current time step.
  //   q: Configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  double terminalCost(Robot& robot, const double t, const Eigen::VectorXd& q, 
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
                             const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  int dimq_, dimv_;
  Eigen::VectorXd lq_, lv_;
  Eigen::VectorXd q_tmp_, v_tmp_;
};

} // namespace idocp


#endif // IDOCP_TERMINAL_PARNMPC_HPP_