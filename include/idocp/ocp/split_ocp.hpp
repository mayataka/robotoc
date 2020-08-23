#ifndef IDOCP_SPLIT_OCP_HPP_
#define IDOCP_SPLIT_OCP_HPP_

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
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
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
  SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
           const std::shared_ptr<Constraints>& constraints);

  // Default constructor.
  SplitOCP();

  // Destructor.
  ~SplitOCP();

  // Use default copy constructor.
  SplitOCP(const SplitOCP&) = default;

  // Use default copy assign operator.
  SplitOCP& operator=(const SplitOCP&) = default;

  // Use default move constructor.
  SplitOCP(SplitOCP&&) noexcept = default;

  // Use default move assign operator.
  SplitOCP& operator=(SplitOCP&&) noexcept = default;
 
  // Check whether the solution s is feasible under inequality constraints.
  bool isFeasible(const Robot& robot, const SplitSolution& s);

  // Initialize the constraints, i.e., set slack and dual variables.
  void initConstraints(const Robot& robot, const int time_step, 
                       const double dtau, const SplitSolution& s);

  // Linearize the OCP for Newton's method around the current solution.
  // Argments: 
  //   robot: The robot model. The contact status of the current time step 
  //      is included in this model.
  //   t: Time of the current time step.
  //   dtau: Discretization length of the OCP.
  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const SplitSolution& s, const SplitSolution& s_next);

  void backwardRiccatiRecursion(const double dtau, 
                                const RiccatiFactorization& riccati_next,
                                RiccatiFactorization& riccati);

  void forwardRiccatiRecursion(const double dtau, SplitDirection& d,   
                               SplitDirection& d_next);

  void computeCondensedDirection(const double dtau, SplitDirection& d);
 
  double maxPrimalStepSize();

  double maxDualStepSize();

  std::pair<double, double> costAndViolation(Robot& robot, const double t, 
                                             const double dtau, 
                                             const SplitSolution& s);

  std::pair<double, double> costAndViolation(Robot& robot, 
                                             const double step_size, 
                                             const double t, const double dtau, 
                                             const SplitSolution& s, 
                                             const SplitSolution& s_next, 
                                             const SplitDirection& d_next);

  void updateDual(const double step_size);

  void updatePrimal(Robot& robot, const double step_size, const double dtau, 
                    const RiccatiFactorization& riccati, 
                    const SplitDirection& d, SplitSolution& s);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) const;

  double squaredKKTErrorNorm(Robot& robot, const double t, const double dtau, 
                             const SplitSolution& s,
                             const SplitSolution& s_next);

private:
  std::shared_ptr<CostFunction> cost_;
  CostFunctionData cost_data_;
  std::shared_ptr<Constraints> constraints_;
  ConstraintsData constraints_data_;
  KKTResidual kkt_residual_;
  KKTMatrix kkt_matrix_;
  StateEquation state_equation_;
  InverseDynamics inverse_dynamics_;
  RiccatiMatrixFactorizer riccati_factorizer_;
  RiccatiMatrixInverter riccati_inverter_;
  Eigen::MatrixXd Ginv_;
  SplitSolution s_tmp_;
  int dimv_, dimf_, dimc_;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = robot.dim_passive() + robot.dimf();
  }

};

} // namespace idocp


#endif // IDOCP_SPLIT_OCP_HPP_