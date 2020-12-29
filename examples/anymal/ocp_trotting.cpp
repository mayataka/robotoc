#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/trotting_configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"


int main(int argc, char *argv[]) {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string path_to_urdf = "../urdf/anymal.urdf";
  idocp::Robot robot(path_to_urdf, contact_frames);

  const double step_length = 0.15;
  const double t_start = 0.5;
  const double t_period = 0.5;

  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << 0, 0, 0.497, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 10, 10, 10, 10, 10, 10, 
              10, 10, 10,
              10, 10, 10,
              10, 10, 10,
              10, 10, 10;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 1, 1, 1, 1, 1, 1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Zero(robot.dimv()));
  a_weight << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01;

  idocp::TrottingSwingAngles swing_angles;
  swing_angles.front_swing_knee   = 1.2; 
  swing_angles.hip_swing_knee     = 1.2;
  auto config_cost = std::make_shared<idocp::TrottingConfigurationSpaceCost>(robot);
  config_cost->set_ref(t_start, t_period, q_standing, step_length, swing_angles);
  config_cost->set_q_weight(q_weight);
  config_cost->set_qf_weight(q_weight);
  config_cost->set_qi_weight(q_weight);
  config_cost->set_v_weight(v_weight);
  config_cost->set_vf_weight(v_weight);
  config_cost->set_vi_weight(v_weight);
  config_cost->set_a_weight(a_weight);
  config_cost->set_dvi_weight(a_weight);
  cost->push_back(config_cost);

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
  contact_cost->set_fi_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(contact_cost);

  auto constraints          = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);

  const double T = 6.05;
  const int N = 240;
  const int max_num_impulse_phase = 11;

  const int num_proc = 4;
  const double t = 0;
  idocp::OCPSolver ocp(robot, cost, constraints, T, N, max_num_impulse_phase, num_proc);

  robot.updateFrameKinematics(q_standing);
  std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  robot.getContactPoints(contact_points);
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  contact_status_initial.setContactPoints(contact_points);
  ocp.setContactStatusUniformly(contact_status_initial);

  auto contact_status_even = robot.createContactStatus();
  contact_status_even.activateContacts({1, 2});
  contact_status_even.setContactPoints(contact_points);
  ocp.pushBackContactStatus(contact_status_even, t_start, t);

  auto contact_status_odd = robot.createContactStatus();
  contact_points[0].coeffRef(0) += 0.5 * step_length;
  contact_points[3].coeffRef(0) += 0.5 * step_length;
  contact_status_odd.activateContacts({0, 3});
  contact_status_odd.setContactPoints(contact_points);
  ocp.pushBackContactStatus(contact_status_odd, t_start+t_period, t);

  for (int i=2; i<=max_num_impulse_phase; ++i) {
    if (i % 2 == 0) {
      contact_points[1].coeffRef(0) += step_length;
      contact_points[2].coeffRef(0) += step_length;
      contact_status_even.setContactPoints(contact_points);
      ocp.pushBackContactStatus(contact_status_even, t_start+i*t_period, t);
    }
    else {
      contact_points[0].coeffRef(0) += step_length;
      contact_points[3].coeffRef(0) += step_length;
      contact_status_odd.setContactPoints(contact_points);
      ocp.pushBackContactStatus(contact_status_odd, t_start+i*t_period, t);
    }
  }

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, // LF
       -0.1, -0.7,  1.0, // LH
        0.1,  0.7, -1.0, // RF
        0.1, -0.7,  1.0; // RH
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp.setStateTrajectory(t, q, v);
  ocp.computeKKTResidual(t, q, v);
  std::cout << "Initial KKT error = " << ocp.KKTError() << std::endl;
  for (int i=0; i<50; ++i) {
    ocp.updateSolution(t, q, v);
  }
  ocp.computeKKTResidual(t, q, v);
  std::cout << "KKT error after iterations = " << ocp.KKTError() << std::endl;

  const std::string path_to_result = "../sim_result/anymal_traj.dat";
  ocp.saveSolution(path_to_result, "q");

  return 0;
}