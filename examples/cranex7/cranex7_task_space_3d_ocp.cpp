#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/unocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/utils/joint_constraints_factory.hpp"


int main() {
  // Create robot
  const std::string path_to_urdf = "../urdf/crane_x7.urdf";
  idocp::Robot robot(path_to_urdf);

  // Create cost function
  auto cost = std::make_shared<idocp::CostFunction>();
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  const int end_effector_frame_id = 26;
  auto task_3d_cost = std::make_shared<idocp::TaskSpace3DCost>(robot, end_effector_frame_id);
  task_3d_cost->set_q_3d_weight(Eigen::Vector3d::Constant(1000));
  task_3d_cost->set_qf_3d_weight(Eigen::Vector3d::Constant(1000));
  Eigen::Vector3d q_ref;
  q_ref << 0, 0, 0.3;
  task_3d_cost->set_q_3d_ref(q_ref);
  cost->push_back(config_cost);
  cost->push_back(task_3d_cost);

  // Create constraints
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();

  // Create OCP solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::UnOCPSolver ocp_solver(robot, cost, constraints, T, N, num_proc);

  ocp_solver.setStateTrajectory(t, q, v);
  ocp_solver.computeKKTResidual(t, q, v);
  std::cout << "Initial KKT error = " << ocp_solver.KKTError() << std::endl;
  const int itr = 30;
  for (int i=0; i<itr; ++i) {
    ocp_solver.updateSolution(t, q, v);
    ocp_solver.computeKKTResidual(t, q, v);
    std::cout << "KKT error after " << i << " iterations = " << ocp_solver.KKTError() << std::endl;
  }
  ocp_solver.printSolution("end-effector", {end_effector_frame_id});
  return 0;
}