#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include "../ocp_benchmarker.hpp"

namespace ocpbenchmark {
namespace iiwa14 {

void BenchmarkWithoutContacts() {
  srand((unsigned int) time(0));
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(joint_cost);
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
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for iiwa 14 without contacts",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.testConvergence(t, q, v, 20);
  ocp_benchmarker.testCPUTime(t, q, v);
  idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for iiwa 14 without contacts",
                                                            robot, cost, constraints, T, N, num_proc);
  parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  parnmpc_benchmarker.testConvergence(t, q, v, 20);
  parnmpc_benchmarker.testCPUTime(t, q, v);
}


void BenchmarkWithContacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {18};
  const double baumgarte_weight_on_velocity = 50;
  const double baumgarte_weight_on_position = 2500;
  const std::string urdf_file_name = "../urdf/iiwa14.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  auto cost = std::make_shared<idocp::CostFunction>();
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->push_back(joint_cost);
  auto constraints = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  const double T = 1;
  const int N = 20;
  const int num_proc = 4;
  const double t = 0;
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  robot.setContactStatus(std::vector<bool>({true}));
  robot.updateKinematics(q, v, Eigen::VectorXd::Zero(robot.dimv()));
  robot.setContactPointsByCurrentKinematics();
  idocp::OCPBenchmarker<idocp::OCP> ocp_benchmarker("OCP for iiwa 14 with contacts",
                                                    robot, cost, constraints, T, N, num_proc);
  ocp_benchmarker.setInitialGuessSolution(t, q, v);
  ocp_benchmarker.setContactStatus(std::vector<bool>({true}));
  ocp_benchmarker.testConvergence(t, q, v, 20);
  ocp_benchmarker.testCPUTime(t, q, v);
  idocp::OCPBenchmarker<idocp::ParNMPC> parnmpc_benchmarker("ParNMPC for iiwa 14 with contacts",
                                                            robot, cost, constraints, T, N, num_proc);
  parnmpc_benchmarker.setInitialGuessSolution(t, q, v);
  parnmpc_benchmarker.setContactStatus(std::vector<bool>({true}));
  parnmpc_benchmarker.testConvergence(t, q, v, 20);
  parnmpc_benchmarker.testCPUTime(t, q, v);
}

} // namespace iiwa14
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::iiwa14::BenchmarkWithoutContacts();
  ocpbenchmark::iiwa14::BenchmarkWithContacts();
  return 0;
}
