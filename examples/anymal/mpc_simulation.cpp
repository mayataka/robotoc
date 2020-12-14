#include <iostream>
#include <string>
#include <memory>
#include <chrono>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
// #include "idocp/ocp/parnmpc.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/time_varying_configuration_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"
#include "idocp/constraints/contact_normal_force.hpp"

#include "idocp/utils/joint_constraints_factory.hpp"
#include "idocp/utils/quadruped_simulator.hpp"

namespace mpcsimulation {
namespace anymal {

void SimulateWithContactsByOCP(
    const std::string& path_to_raisim_activation_key) {
  srand((unsigned int) time(0));
  const std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string urdf_file_name = "../anymal/anymal.urdf";
  idocp::Robot robot(urdf_file_name, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_ref(robot.dimq());
  q_ref << 0, 0, 0.4792, 0, 0, 0, 1, 
           -0.1,  0.7, -1.0, 
           -0.1, -0.7,  1.0, 
            0.1,  0.7, -1.0, 
            0.1, -0.7,  1.0;
  robot.normalizeConfiguration(q_ref);
  Eigen::VectorXd v_ref(robot.dimv());
  v_ref << 0, 0, 0, 0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0;
  auto configuration_cost = std::make_shared<idocp::TimeVaryingConfigurationCost>(robot);
  configuration_cost->set_ref(0, q_ref, v_ref);
  Eigen::VectorXd q_weight(robot.dimv());
  // q_weight << 0, 0, 0, 0, 0, 0, 
  // q_weight << 1000, 1000, 1000, 1000, 1000, 1000, 
              // 1000, 1000, 1000, 
              // 1000, 1000, 1000, 
              // 1000, 1000, 1000, 
              // 1000, 1000, 1000;
  // configuration_cost->set_q_weight(q_weight);
  // configuration_cost->set_qf_weight(q_weight);
  configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  configuration_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  configuration_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  configuration_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  cost->push_back(configuration_cost);

  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw;
    // fw << 1, 1, 0.001;
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
  }
  contact_cost->set_f_weight(f_weight);
  cost->push_back(contact_cost);
  auto impulse_configuration_cost = std::make_shared<idocp::ImpulseTimeVaryingConfigurationCost>(robot);
  impulse_configuration_cost->set_ref(0, q_ref, v_ref);
  // impulse_configuration_cost->set_q_weight(q_weight);
  impulse_configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  impulse_configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  impulse_configuration_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  cost->push_back(impulse_configuration_cost);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
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
  const double T = 1;
  const int N = 20;
  const int max_num_impulse_phase = 4;
  const int num_proc = 4;
  const double t = 0;
  Eigen::VectorXd q(robot.dimq());
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, 
       -0.1, -0.7,  1.0, 
        0.1,  0.7, -1.0, 
        0.1, -0.7,  1.0;
  Eigen::VectorXd v(robot.dimv());
  v << 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 
       0, 0, 0, 
       0, 0, 0, 
       0, 0, 0;
  auto contact_status = robot.createContactStatus();
  robot.updateFrameKinematics(q);
  robot.setContactPoints(contact_status);
  contact_status.activateContacts({0, 1, 2, 3});

  idocp::MPC<idocp::OCP> mpc(robot, cost, constraints, T, N, max_num_impulse_phase, num_proc);
  mpc.getSolverHandle()->setContactStatusUniformly(contact_status);
  mpc.initializeSolution(t, q, v, 100);
  mpc.getSolverHandle()->printSolution();

  const std::string urdf_for_raisim_file_name = "/home/sotaro/src/idocp/examples/anymal/anymal/anymal_for_raisim.urdf";
  idocp::QuadrupedSimulator simulator(path_to_raisim_activation_key, 
                                      urdf_for_raisim_file_name, 
                                      "../sim_result", "forward");
  constexpr bool visualization = true;
  constexpr bool video_recording = true;
  simulator.run(robot, mpc, 10, 0.0025, 0, q, v, visualization, video_recording);
}

} // namespace anymal 
} // namespace ocpbenchmark


int main(int argc, char *argv[]) {
  mpcsimulation::anymal::SimulateWithContactsByOCP(argv[1]);
  return 0;
}
