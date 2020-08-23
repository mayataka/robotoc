#ifndef IDOCP_TERMINAL_OCP_HPP_
#define IDOCP_TERMINAL_OCP_HPP_

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
#include "riccati_factorization.hpp"


namespace idocp {

class TerminalOCP {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Sets the robot, cost function, and constraints.
  // Argments:
  //    robot: The robot model that has been already initialized.
  //    cost: The pointer to the cost function.
  //    constraints: The pointer to the constraints.
  TerminalOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints);

  // Default constructor.
  TerminalOCP();
  
  // Destructor.
  ~TerminalOCP();

  // Use default copy constructor.
  TerminalOCP(const TerminalOCP&) = default;

  // Use default copy assign operator.
  TerminalOCP& operator=(const TerminalOCP&) = default;

  // Use default move constructor.
  TerminalOCP(TerminalOCP&&) noexcept = default;

  // Use default move assign operator.
  TerminalOCP& operator=(TerminalOCP&&) noexcept = default;

  // Check whether the solution s is feasible under inequality constraints.
  bool isFeasible(const Robot& robot, const SplitSolution& s);

  // Initialize the constraints, i.e., set slack and dual variables.
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  // Linearize the OCP for Newton's method around the current solution at the 
  // last time step of the horizon.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   dtau: Discretization interval of the horizon.
  void linearizeOCP(Robot& robot, const double t, const SplitSolution& s, 
                    RiccatiFactorization& riccati);

  // Computes the direction of the condensed veriables.
  // Argments: 
  //   dtau: Discretization interval of the horizon.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  void computeCondensedDirection(Robot& robot, const double dtau, 
                                 SplitDirection& d);

  // Returns the maximum step size of the primal variables of the inequality 
  // constraints.
  double maxPrimalStepSize();

  // Returns the maximum step size of the dual variables of the inequality 
  // constraints.
  double maxDualStepSize();

  // Returns the terminal cost.
  // Argments: 
  //   t: Time of the current time step.
  double terminalCost(Robot& robot, const double t, const SplitSolution& s);

  // Returns the terminal cost with step_size.
  // Argments: 
  //   robot: The robot model. 
  //   step_size: The step size.
  //   t: Time of the current time step.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  double terminalCost(Robot& robot, const double step_size, const double t, 
                      const SplitSolution& s, const SplitDirection& d);

  // Updates the dual variables of the inequality constraints.
  // Argments: 
  //   step_size: The step size.
  void updateDual(const double step_size);

  // Updates the primal variables.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   step_size: The step size.
  //   dq: Direction of the configuration. Size must be dimv.
  //   dv: Direction of the generalized velocity. Size must be dimv.
  void updatePrimal(Robot& robot, const double step_size, 
                    const RiccatiFactorization& riccati,
                    const SplitDirection& d, SplitSolution& s) const;

  // Returns the squared KKT error norm.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  double squaredKKTErrorNorm(Robot& robot, const double t, 
                             const SplitSolution& s);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  SplitSolution s_tmp_;

};

} // namespace idocp


#endif // IDOCPT_TERMINAL_OCP_HPP_