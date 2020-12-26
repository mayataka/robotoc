#include <iostream>
#include <string>
#include <memory>
#include <deque>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
// #include "idocp/ocp/parnmpc_solver.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/cost/time_varying_configuration_cost.hpp"
#include "idocp/cost/trotting_configuration_cost.hpp"
#include "idocp/cost/trotting_foot_step_cost.hpp"
#include "idocp/cost/impulse_trotting_foot_step_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/cost/foot_step_trotting_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"
#include "idocp/constraints/contact_normal_force.hpp"
#include "idocp/constraints/contact_distance.hpp"

#include "idocp/utils/trajectory_viewer.hpp"



int main(int argc, char *argv[]) {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string path_to_urdf = "../anymal/anymal.urdf";
  idocp::Robot robot(path_to_urdf, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_ref(Eigen::VectorXd::Zero(robot.dimq()));
  q_ref <<  0, 0, 0.4792, 0, 0, 0, 1, 
           -0.1,  0.7, -1.0, 
           -0.1, -0.7,  1.0, 
            0.1,  0.7, -1.0, 
            0.1, -0.7,  1.0;
  Eigen::VectorXd v_ref(Eigen::VectorXd::Zero(robot.dimv()));
  v_ref <<  0.075, 0, 0, 0, 0, 0, 
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, 0, 0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  // q_weight << 0, 0, 0, 0, 0, 0, 
  q_weight << 10, 10, 10, 10, 10, 10, 
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 1, 1, 1, 1, 1, 1, 
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Zero(robot.dimv()));
  a_weight << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_a_weight(a_weight);
  joint_cost->set_qf_weight(q_weight);
  joint_cost->set_vf_weight(v_weight);
  cost->push_back(joint_cost);

  auto LF_foot_cost = std::make_shared<idocp::TrottingFootStepCost>(robot, 14);
  auto LH_foot_cost = std::make_shared<idocp::TrottingFootStepCost>(robot, 24);
  auto RF_foot_cost = std::make_shared<idocp::TrottingFootStepCost>(robot, 34);
  auto RH_foot_cost = std::make_shared<idocp::TrottingFootStepCost>(robot, 44);
  robot.updateFrameKinematics(q_ref);
  Eigen::Vector3d LF_ref = robot.framePosition(14);
  Eigen::Vector3d LH_ref = robot.framePosition(24);
  Eigen::Vector3d RF_ref = robot.framePosition(34);
  Eigen::Vector3d RH_ref = robot.framePosition(44);
  // LF_ref.coeffRef(0) += 0.5 * 0.15; RH_ref.coeffRef(0) += 0.5 * 0.15;
  // LF_foot_cost->set_q_3d_ref(LF_ref);
  // LH_foot_cost->set_q_3d_ref(LH_ref);
  // RF_foot_cost->set_q_3d_ref(RF_ref);
  // RH_foot_cost->set_q_3d_ref(RH_ref);
  const double step_length = 0.15;
  const double step_height = 0.0;
  LF_foot_cost->set_q_3d_ref(LF_ref, step_length, step_height);
  LH_foot_cost->set_q_3d_ref(LH_ref, step_length, step_height);
  RF_foot_cost->set_q_3d_ref(RF_ref, step_length, step_height);
  RH_foot_cost->set_q_3d_ref(RH_ref, step_length, step_height);
  const double t_start = 0.5;
  const double t_period = 0.5;
  LF_foot_cost->set_period(t_start,          2*t_period, true);
  LH_foot_cost->set_period(t_start+t_period, 2*t_period, false);
  RF_foot_cost->set_period(t_start+t_period, 2*t_period, false);
  RH_foot_cost->set_period(t_start,          2*t_period, true);
  const Eigen::Vector3d q_3d_weight = Eigen::Vector3d::Constant(1000);
  LF_foot_cost->set_q_3d_weight(q_3d_weight); LF_foot_cost->set_qf_3d_weight(q_3d_weight);
  LH_foot_cost->set_q_3d_weight(q_3d_weight); LH_foot_cost->set_qf_3d_weight(q_3d_weight);
  RF_foot_cost->set_q_3d_weight(q_3d_weight); RF_foot_cost->set_qf_3d_weight(q_3d_weight);
  RH_foot_cost->set_q_3d_weight(q_3d_weight); RH_foot_cost->set_qf_3d_weight(q_3d_weight);
  // cost->push_back(LF_foot_cost);
  // cost->push_back(LH_foot_cost);
  // cost->push_back(RF_foot_cost);
  // cost->push_back(RH_foot_cost);

  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    fr << 0, 0, 70;
    // fr << 0, 0, 0;
    f_ref.push_back(fr);
  }
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(contact_cost);
  auto impulse_joint_cost = std::make_shared<idocp::JointSpaceImpulseCost>(robot);
  impulse_joint_cost->set_q_ref(q_ref);
  impulse_joint_cost->set_v_ref(v_ref);
  impulse_joint_cost->set_q_weight(q_weight);
  impulse_joint_cost->set_v_weight(v_weight);
  impulse_joint_cost->set_dv_weight(a_weight);
  cost->push_back(impulse_joint_cost);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  impulse_force_cost->set_f_ref(f_ref);
  cost->push_back(impulse_force_cost);

  auto impulse_LF_foot_cost = std::make_shared<idocp::ImpulseTrottingFootStepCost>(robot, 14);
  auto impulse_LH_foot_cost = std::make_shared<idocp::ImpulseTrottingFootStepCost>(robot, 24);
  auto impulse_RF_foot_cost = std::make_shared<idocp::ImpulseTrottingFootStepCost>(robot, 34);
  auto impulse_RH_foot_cost = std::make_shared<idocp::ImpulseTrottingFootStepCost>(robot, 44);
  impulse_LF_foot_cost->set_q_3d_ref(LF_ref, step_length, step_height);
  impulse_LH_foot_cost->set_q_3d_ref(LH_ref, step_length, step_height);
  impulse_RF_foot_cost->set_q_3d_ref(RF_ref, step_length, step_height);
  impulse_RH_foot_cost->set_q_3d_ref(RH_ref, step_length, step_height);
  impulse_LF_foot_cost->set_period(t_start,          2*t_period, true);
  impulse_LH_foot_cost->set_period(t_start+t_period, 2*t_period, false);
  impulse_RF_foot_cost->set_period(t_start+t_period, 2*t_period, false);
  impulse_RH_foot_cost->set_period(t_start,          2*t_period, true);
  impulse_LF_foot_cost->set_q_3d_weight(q_3d_weight); impulse_LF_foot_cost->set_qf_3d_weight(q_3d_weight);
  impulse_LH_foot_cost->set_q_3d_weight(q_3d_weight); impulse_LH_foot_cost->set_qf_3d_weight(q_3d_weight);
  impulse_RF_foot_cost->set_q_3d_weight(q_3d_weight); impulse_RF_foot_cost->set_qf_3d_weight(q_3d_weight);
  impulse_RH_foot_cost->set_q_3d_weight(q_3d_weight); impulse_RH_foot_cost->set_qf_3d_weight(q_3d_weight);
  cost->push_back(impulse_LF_foot_cost);
  cost->push_back(impulse_LH_foot_cost);
  cost->push_back(impulse_RF_foot_cost);
  cost->push_back(impulse_RH_foot_cost);

  auto constraints          = std::make_shared<idocp::Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  auto contact_normal_force = std::make_shared<idocp::ContactNormalForce>(robot);
  auto contact_distance     = std::make_shared<idocp::ContactDistance>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  // constraints->push_back(contact_normal_force);
  // constraints->push_back(contact_distance);
  const double T = 1.0;
  // const double T = 0.75;
  // const double T = 2.5;
  const int N = 40;
  // const double T = 1.5;
  // const int N = 60;
  const int max_num_impulse_phase = 4;
  const int num_proc = 4;
  const double t = 0;
  idocp::OCPSolver ocp(robot, cost, constraints, T, N, max_num_impulse_phase, num_proc);

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, // LF
       -0.1, -0.7,  1.0, // LH
        0.1,  0.7, -1.0, // RF
        0.1, -0.7,  1.0; // RH
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  robot.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  robot.getContactPoints(contact_points);
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  contact_status_initial.setContactPoints(contact_points);
  ocp.setContactStatusUniformly(contact_status_initial);

  // contact_points[0].coeffRef(0) += 0.5 * step_length;
  // contact_points[3].coeffRef(0) += 0.5 * step_length;
  auto contact_status_even = robot.createContactStatus();
  // contact_status_even.activateContacts({0, 3});
  contact_status_even.activateContacts({1, 2});
  contact_status_even.setContactPoints(contact_points);
  ocp.pushBackContactStatus(contact_status_even, 0.5, t);


  // auto contact_status_initial = robot.createContactStatus();
  // contact_status_initial.activateContacts({0, 1, 2, 3});
  // robot.updateFrameKinematics(q);
  // std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  // robot.getContactPoints(contact_points);
  // contact_status_initial.setContactPoints(contact_points);
  // ocp.setContactStatusUniformly(contact_status_initial);

  // auto contact_status_even = robot.createContactStatus();
  // contact_status_even.activateContacts({1, 2});
  // contact_status_even.setContactPoints(contact_points);
  // ocp.pushBackContactStatus(contact_status_even, t_start, t);

  // auto contact_status_odd = robot.createContactStatus();
  // contact_points[0].coeffRef(0) += 0.5 * step_length;
  // contact_points[3].coeffRef(0) += 0.5 * step_length;
  // contact_status_odd.activateContacts({0, 3});
  // contact_status_odd.setContactPoints(contact_points);
  // ocp.pushBackContactStatus(contact_status_odd, t_start+t_period, t);

  ocp.setStateTrajectory(t, q, v);
  ocp.computeKKTResidual(t, q, v);
  std::cout << "Initial KKT error = " << ocp.KKTError() << std::endl;
  for (int i=0; i<10; ++i) {
    ocp.updateSolution(t, q, v);
  }
  ocp.computeKKTResidual(t, q, v);
  std::cout << "KKT error after iterations = " << ocp.KKTError() << std::endl;
  ocp.printSolution("end-effector", {14, 24, 34, 44});

  // const double dt = 0.1;
  // for (int i=0; i<100; ++i) {
  //   std::cout << "LF_ref[" << t+i*dt << "] = " << LF_foot_cost->q_3d_ref(t+i*dt).transpose() << std::endl;
  //   std::cout << "LH_ref[" << t+i*dt << "] = " << LH_foot_cost->q_3d_ref(t+i*dt).transpose() << std::endl;
  //   std::cout << "RF_ref[" << t+i*dt << "] = " << RF_foot_cost->q_3d_ref(t+i*dt).transpose() << std::endl;
  //   std::cout << "RH_ref[" << t+i*dt << "] = " << RH_foot_cost->q_3d_ref(t+i*dt).transpose() << std::endl;
  // }

  const std::string path_to_result = "../sim_result/anymal_q.dat";
  ocp.saveSolution(path_to_result, "q");
  const std::string path_to_raisim_activation_key = argv[1];
  const std::string path_to_urdf_for_raisim = "../anymal/anymal_for_raisim.urdf";
  idocp::TrajectoryViewer viewer(path_to_raisim_activation_key, path_to_urdf_for_raisim);
  viewer.display(path_to_result, T/N);
  return 0;
}