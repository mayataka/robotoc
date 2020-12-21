#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
// #include "idocp/ocp/parnmpc_solver.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/time_varying_configuration_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/ocp_benchmarker.hpp"


int main () {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 0.4792, 0, 0, 0, 1, 
           -0.1,  0.7, -1.0, 
           -0.1, -0.7,  1.0, 
            0.1,  0.7, -1.0, 
            0.1, -0.7,  1.0;
  Eigen::VectorXd v_ref(robot.dimv());
  v_ref << 0, 0, 0, 0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0;
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimu(), 0.001));
  cost->push_back(joint_cost);
  // auto configuration_cost = std::make_shared<idocp::TimeVaryingConfigurationCost>(robot);
  // configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // configuration_cost->set_ref(0, q_ref, v_ref);
  // configuration_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // configuration_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // configuration_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  // cost->push_back(configuration_cost);
  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    fr << 0, 0, 70;
    f_ref.push_back(fr);
  }
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(contact_cost);
  auto impulse_joint_cost = std::make_shared<idocp::JointSpaceImpulseCost>(robot);
  impulse_joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  impulse_joint_cost->set_q_ref(q_ref);
  impulse_joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  impulse_joint_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(impulse_joint_cost);
  // auto impulse_configuration_cost = std::make_shared<idocp::ImpulseTimeVaryingConfigurationCost>(robot);
  // impulse_configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // impulse_configuration_cost->set_ref(0, q_ref, v_ref);
  // impulse_configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // impulse_configuration_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  // cost->push_back(impulse_configuration_cost);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  impulse_force_cost->set_f_ref(f_ref);
  cost->push_back(impulse_force_cost);
  auto constraints = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  const double T = 0.5;
  const int N = 20;
  const int max_num_impulse_phase = 5;
  const int num_proc = 1;
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, 
       -0.1, -0.7,  1.0, 
        0.1,  0.7, -1.0, 
        0.1, -0.7,  1.0;
  Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  v << 0, 0, 0, 0, 0, 0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0;
  idocp::OCPBenchmarker<idocp::OCPSolver> ocp_benchmarker("OCP for anymal with contacts",
                                                          robot, cost, constraints, T, N, 
                                                          max_num_impulse_phase, num_proc);
  auto contact_status = robot.createContactStatus();
  // contact_status.activateContacts({0, 3});
  contact_status.activateContacts({0, 1, 2, 3});
  robot.updateFrameKinematics(q);
  robot.setContactPoints(contact_status);
  ocp_benchmarker.getSolverHandle()->setContactStatusUniformly(contact_status);
  auto contact_status_next = robot.createContactStatus();
  // contact_status_next.activateContacts({0, 1, 2, 3});
  contact_status_next.activateContacts({1, 2});
  robot.updateFrameKinematics(q);
  robot.setContactPoints(contact_status_next);
  const double switching_time = 0.18;
  ocp_benchmarker.getSolverHandle()->pushBackContactStatus(contact_status_next, switching_time, t);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 50, false);
  ocp_benchmarker.testConvergence(t, q, v, 10, false);
  const double t_next = t + 0.001;
  // ocp_benchmarker.getSolverHandle()->warmStartSolution(t, t_next);
  ocp_benchmarker.testConvergence(t_next, q, v, 10, false);
  // ocp_benchmarker.testCPUTime(t, q, v, 1000);
  // ocp_benchmarker.getSolverHandle()->printSolution("q");
  // ocp_benchmarker.getSolverHandle()->printSolution("v");
  // idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for anymal with contacts",
  //                                                           robot, cost, constraints, T, N, num_proc);
  // parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  // parnmpc_benchmarker.getSolverHandle()->activateContacts({0, 1, 2, 3}, 0, N);
  // parnmpc_benchmarker.testConvergence(t, q, v, 20, false);
  // parnmpc_benchmarker.testCPUTime(t, q, v);

  return 0;
}
