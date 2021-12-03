#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/unconstr/unconstr_ocp.hpp"
#include "robotoc/solver/unconstr_ocp_solver.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/utils/joint_constraints_factory.hpp"
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
  robotoc::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();
  constraints->setBarrier(1.0e-03);

  // Create the OCP solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  const int nthreads = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Constant(robot.dimq(), 2);
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  robotoc::UnconstrOCP ocp(robot, cost, constraints, T, N);
  robotoc::UnconstrOCPSolver ocp_solver(ocp, nthreads);

  // Solves the OCP.
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  const int num_iteration = 30;
  const bool line_search = false;
  robotoc::benchmark::convergence(ocp_solver, t, q, v, num_iteration, line_search);
  const int num_iteration_CPU = 10000;
  robotoc::benchmark::CPUTime(ocp_solver, t, q, v, num_iteration_CPU, line_search);

  return 0;
}