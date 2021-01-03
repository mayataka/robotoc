#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/unparnmpc_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/task_space_6d_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/ocp_benchmarker.hpp"


int main() {
  // Creates a robot.
  const std::string path_to_urdf = "../urdf/crane_x7.urdf";
  idocp::Robot robot(path_to_urdf);
  
  // Creates a cost function.
  auto cost = std::make_shared<idocp::CostFunction>();
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  config_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  const int end_effector_frame_id = 26;
  auto task_6d_cost = std::make_shared<idocp::TaskSpace6DCost>(robot, end_effector_frame_id);
  task_6d_cost->set_q_6d_weight(Eigen::VectorXd::Constant(3, 1000), Eigen::VectorXd::Constant(3, 1000));
  task_6d_cost->set_qf_6d_weight(Eigen::VectorXd::Constant(3, 1000), Eigen::VectorXd::Constant(3, 1000));
  Eigen::Vector3d q_position_ref;
  q_position_ref << 0, 0, 0.3;
  const Eigen::Matrix3d q_rotation_ref = Eigen::Matrix3d::Identity();
  task_6d_cost->set_q_6d_ref(q_position_ref, q_rotation_ref);
  cost->push_back(config_cost);
  cost->push_back(task_6d_cost);

  // Creates joint constraints.
  idocp::JointConstraintsFactory constraints_factory(robot);
  auto constraints = constraints_factory.create();

  // Creates the ParNMPC solver for unconstrained rigid-body systems.
  const double T = 1;
  const int N = 20;
  const int nthreads = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  idocp::UnParNMPCSolver parnmpc_solver(robot, cost, constraints, T, N, nthreads);

  // Solves the OCP.
  parnmpc_solver.setStateTrajectory(t, q, v);
  parnmpc_solver.initBackwardCorrection(t);
  const int num_iteration = 30;
  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(parnmpc_solver, t, q, v, num_iteration, line_search);
  parnmpc_solver.printSolution("end-effector", {end_effector_frame_id});

  return 0;
}