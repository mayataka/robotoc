#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/unconstr/unconstr_parnmpc.hpp"
#include "robotoc/solver/unconstr_parnmpc_solver.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/utils/ocp_benchmarker.hpp"


int main() {
  // Create a robot.
  const std::string path_to_urdf = "../iiwa_description/urdf/iiwa14.urdf";
  robotoc::Robot robot(path_to_urdf);

  // Create a cost function.
  robot.setJointEffortLimit(Eigen::VectorXd::Constant(robot.dimu(), 200));
  auto cost = std::make_shared<robotoc::CostFunction>();
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_ref(Eigen::VectorXd::Constant(robot.dimv(), -5));
  config_cost->set_v_ref(Eigen::VectorXd::Constant(robot.dimv(), -9));
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  config_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.0));
  cost->push_back(config_cost);

  // Create joint constraints.
  const double barrier = 1.0e-03;
  const double fraction_to_boundary_rule = 0.995;
  auto constraints = std::make_shared<robotoc::Constraints>(barrier, fraction_to_boundary_rule);
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);

  // Create the ParNMPC solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  robotoc::UnconstrParNMPC parnmpc(robot, cost, constraints, T, N);
  auto solver_options = robotoc::SolverOptions::defaultOptions();
  const int nthreads = 8; // Please set nthreads by the number of the processors of your PC to enjoy ParNMPC!
  robotoc::UnconstrParNMPCSolver parnmpc_solver(parnmpc, solver_options, nthreads);

  // Initial time and initial state
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Constant(robot.dimq(), 2);
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());

  // Solves the OCP.
  parnmpc_solver.setSolution("q", q);
  parnmpc_solver.setSolution("v", v);
  parnmpc_solver.initConstraints();
  parnmpc_solver.initBackwardCorrection(t);
  std::cout << "Initial KKT error: " << parnmpc_solver.KKTError(t, q, v) << std::endl;
  parnmpc_solver.solve(t, q, v);
  std::cout << "KKT error after convergence: " << parnmpc_solver.KKTError(t, q, v) << std::endl;
  std::cout << parnmpc_solver.getSolverStatistics() << std::endl;

  // Measures CPU timing
  const int num_iteration_CPU = 10000;
  robotoc::benchmark::CPUTime(parnmpc_solver, t, q, v, num_iteration_CPU);

  return 0;
}