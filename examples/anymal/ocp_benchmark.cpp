#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"


namespace ocpbenchmark {
namespace anymal {

void CPUTime_with_contacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baumgarte_weight_on_velocity = 10;
  const double baumgarte_weight_on_position = 10;
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 4.8, 0, 0, 0, 1, 
           0.0315, 0.4, -0.8, 
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
  const int N = 50;
  const int num_proc = 4;
  idocp::ParNMPC parnmpc(robot, cost, constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 4.8, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  parnmpc.setStateTrajectory(q, v);
  std::vector<bool> contact_status = {true, true, true, true};
  std::vector<std::vector<bool>> contact_sequence = {N, contact_status};
  parnmpc.setContactSequence(contact_sequence);
  const int num_iteration = 100;
  std::chrono::system_clock::time_point start_clock, end_clock;
  start_clock = std::chrono::system_clock::now();
  for (int i=0; i<num_iteration; ++i) {
    parnmpc.updateSolution(t, q, v, false);
  }
  end_clock = std::chrono::system_clock::now();
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "model: anymal" << std::endl;
  std::cout << "dimq = " << robot.dimq() << std::endl;
  std::cout << "dimv = " << robot.dimv() << std::endl;
  std::cout << "max_dimf = " << robot.max_dimf() << std::endl;
  std::cout << "N = " << N << std::endl;
  std::cout << "T = " << T << std::endl;
  std::cout << "number of threads = " << num_proc << std::endl;
  std::cout << "total CPU time: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() << "[ms]" << std::endl;
  std::cout << "CPU time per update: " << 1e-03 * std::chrono::duration_cast<std::chrono::microseconds>(end_clock-start_clock).count() / num_iteration << "[ms]" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}


void KKTError_with_contacts() {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const double baumgarte_weight_on_velocity = 10;
  const double baumgarte_weight_on_position = 100;
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames, 
                     baumgarte_weight_on_velocity, 
                     baumgarte_weight_on_position);
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 4.8, 0, 0, 0, 1, 
           0.0315, 0.4, -0.8, 
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
  const int N = 50;
  const int num_proc = 4;
  idocp::ParNMPC parnmpc(robot, cost, constraints, T, N, num_proc);
  const double t = 0;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  q << 0, 0, 4.8, 0, 0, 0, 1, 
       0.0315, 0.4, -0.8, 
       0.0315, -0.4, 0.8, 
       -0.0315, 0.4, -0.8,
       -0.0315, -0.4, 0.8;
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(robot.dimv());
  parnmpc.setStateTrajectory(q, v);
  std::vector<bool> contact_status = {true, true, true, true};
  std::vector<std::vector<bool>> contact_sequence = {N, contact_status};
  parnmpc.setContactSequence(contact_sequence);
  std::cout << "---------- OCP benchmark ----------" << std::endl;
  std::cout << "model: anymal" << std::endl;
  std::cout << "dimq = " << robot.dimq() << std::endl;
  std::cout << "dimv = " << robot.dimv() << std::endl;
  std::cout << "max_dimf = " << robot.max_dimf() << std::endl;
  std::cout << "q = " << q.transpose() << std::endl;
  std::cout << "v = " << v.transpose() << std::endl;
  std::cout << "Initial KKT error = " << parnmpc.KKTError(t, q, v) << std::endl;
  const int num_iteration = 50;
  for (int i=0; i<num_iteration; ++i) {
    parnmpc.updateSolution(t, q, v, true);
    std::cout << "KKT error at iteration " << i << " = " << parnmpc.KKTError(t, q, v) << std::endl;
  }
  std::cout << "-----------------------------------" << std::endl;
  std::cout << std::endl;
}

} // namespace anymal
} // namespace ocpbenchmark


int main() {
  ocpbenchmark::anymal::CPUTime_with_contacts();
  ocpbenchmark::anymal::KKTError_with_contacts();
  return 0;
}
