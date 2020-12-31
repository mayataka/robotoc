#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/ocp_benchmarker.hpp"


int main() {
  // Create a robot.
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name);

  // Create a cost function.
  robot.setJointEffortLimit(Eigen::VectorXd::Constant(robot.dimu(), 200));
  auto cost = std::make_shared<idocp::CostFunction>();
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
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
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();

  // Create the OCP solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Constant(robot.dimq(), 2);
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, 0, num_proc);

  // Solves the OCP.
  ocp_solver.setStateTrajectory(t, q, v);
  const int num_iteration = 50;
  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, num_iteration, line_search);
  // const int num_iteration_CPU = 10000;
  // idocp::ocpbenchmarker::CPUTime(ocp_solver, t, q, v, num_iteration_CPU, line_search);
  // ocp_solver.printSolution();

  return 0;
}