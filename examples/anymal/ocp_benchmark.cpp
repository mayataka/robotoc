#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
// #include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
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
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/discrete_event.hpp"

namespace ocpbenchmark {
namespace anymal {

void BenchmarkWithContacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames);
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 0.48, 0, 0, 0, 1, 
           0.0315, 0.4, -0.6, 
           0.0315, -0.4, 0.8, 
           -0.0315, 0.4, -0.8,
           -0.0315, -0.4, 0.8;
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimu(), 0.0));
  // joint_cost->set_u_passive_weight(Eigen::VectorXd::Constant(6, 0.001));
  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<contact_frames.size(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_cost->set_f_weight(f_weight);
  auto impulse_joint_cost = std::make_shared<idocp::JointSpaceImpulseCost>(robot);
  impulse_joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  impulse_joint_cost->set_q_ref(q_ref);
  impulse_joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  impulse_joint_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  cost->push_back(joint_cost);
  cost->push_back(contact_cost);
  cost->push_back(impulse_joint_cost);
  cost->push_back(impulse_force_cost);
//   idocp::JointConstraintsFactory constraints_factory(robot);
//   auto constraints = constraints_factory.create();
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
  // constraints->push_back(joint_torques_lower);
  // constraints->push_back(joint_torques_upper);
  const double T = 1;
  const int N = 20;
  const int max_num_impulse_phase = 5;
  const int num_proc = 4;
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 0.48, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  v << 0.0, 0, 0, 0, 0, 0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0, 
       0.0, 0.0, 0.0,
       0.0, 0.0, 0.0;
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for anymal with contacts",
                                                    robot, cost, constraints, T, N, 
                                                    max_num_impulse_phase, num_proc);
  q << 0, 0, 0.5, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  auto contact_status = robot.createContactStatus();
  robot.updateFrameKinematics(q);
  robot.setContactPoints(contact_status);
  // contact_status.activateContacts({0, 1, 2, 3});
  contact_status.activateContacts({1, 2, 3});
  ocp_benchmarker.getSolverHandle()->setContactStatusUniformly(contact_status);
  auto contact_status_next = robot.createContactStatus();
  robot.setContactPoints(contact_status_next);
  // contact_status_next.activateContacts({1, 2, 3});
  contact_status_next.activateContacts({0, 1, 2, 3});
  idocp::DiscreteEvent discrete_event(robot);
  discrete_event.setDiscreteEvent(contact_status, contact_status_next);
  discrete_event.setContactPoints(robot);
  // discrete_event.eventTime = 0.475;
  discrete_event.eventTime = 0.175;
  // discrete_event.eventTime = 0.125;
  // discrete_event.eventTime = 0.025;
  ocp_benchmarker.getSolverHandle()->setDiscreteEvent(discrete_event);
  q << 0, 0, 0.5, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 50, false);
  // ocp_benchmarker.testCPUTime(t, q, v, 10000);
  // idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for anymal with contacts",
  //                                                           robot, cost, constraints, T, N, num_proc);
  // parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  // parnmpc_benchmarker.getSolverHandle()->activateContacts({0, 1, 2, 3}, 0, N);
  // parnmpc_benchmarker.testConvergence(t, q, v, 20, false);
  // parnmpc_benchmarker.testCPUTime(t, q, v);
}

} // namespace anymal
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::anymal::BenchmarkWithContacts();
  return 0;
}
