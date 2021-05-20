#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/solver/unocp_solver.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "idocp/utils/trajectory_viewer.hpp"
#endif 


int main(int argc, char *argv[]) {
  // Create a robot.
  const std::string path_to_urdf = "../iiwa_description/urdf/iiwa14.urdf";
  idocp::Robot robot(path_to_urdf);

  // Change the limits from the default parameters.
  robot.setJointEffortLimit(Eigen::VectorXd::Constant(robot.dimu(), 50));
  robot.setJointVelocityLimit(Eigen::VectorXd::Constant(robot.dimv(), M_PI_2));

  // Create a cost function.
  auto cost = std::make_shared<idocp::CostFunction>();
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_ref(Eigen::VectorXd::Zero(robot.dimq()));
  q_ref << 0, M_PI_2, 0, M_PI_2, 0, M_PI_2, 0;
  config_cost->set_q_ref(q_ref);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(config_cost);

  // Create joint constraints.
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();

  // Create the OCP solver for unconstrained rigid-body systems.
  const double T = 3;
  const int N = 60;
  const int nthreads = 4;
  const double t = 0;
  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << M_PI_2, 0, M_PI_2, 0, M_PI_2, 0, M_PI_2;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::UnOCPSolver ocp_solver(robot, cost, constraints, T, N, nthreads);

  // Solves the OCP.
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  const int num_iteration = 30;
  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, num_iteration, line_search);

#ifdef ENABLE_VIEWER
  idocp::TrajectoryViewer viewer(path_to_urdf);
  const double dt = T/N;
  viewer.display(ocp_solver.getSolution("q"), dt);
#endif 

  return 0;
}