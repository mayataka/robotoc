#include <iostream>
#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/time_varying_configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"

#include "idocp/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "idocp/utils/trajectory_viewer.hpp"
#endif 


int main(int argc, char *argv[]) {
  std::vector<int> contact_frames = {14, 24, 34, 44}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  idocp::Robot robot(path_to_urdf, contact_frames);

  const double stride = 0.4;
  const double additive_stride_hip = 0.2;
  const double t_start = 1.0;

  const double t_front_swing = 0.135;
  const double t_front_hip_swing = 0.05;
  const double t_hip_swing = 0.165;
  const double t_period = t_front_swing + t_front_hip_swing + t_hip_swing;
  const int steps = 10;

  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << -3, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 1, 1, 1, 10, 10, 10, 
              10, 10, 10,
              10, 10, 10,
              10, 10, 10,
              10, 10, 10;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Zero(robot.dimv()));
  a_weight << 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01;

  auto config_cost = std::make_shared<idocp::TimeVaryingConfigurationSpaceCost>(robot);
  Eigen::VectorXd q_ref = q_standing;
  Eigen::VectorXd v_ref(Eigen::VectorXd::Zero(robot.dimv()));
  v_ref(0) = stride / t_period;
  config_cost->set_ref(robot, t_start, t_start+(0.5+steps)*t_period, q_standing, v_ref);
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
    fw << 1e-01, 1e-01, 1.0e-07;
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

  const double T = 7; 
  const int N = 140;
  const int max_num_impulse_phase = (steps+3)*2;

  const int nthreads = 4;
  const double t = 0;
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, max_num_impulse_phase, nthreads);

  robot.updateFrameKinematics(q_standing);
  std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  robot.getContactPoints(contact_points);
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  auto contact_status_front_swing = robot.createContactStatus();
  contact_status_front_swing.activateContacts({1, 3});
  auto contact_status_hip_swing = robot.createContactStatus();
  contact_status_hip_swing.activateContacts({0, 2});
  auto contact_status_front_hip_swing = robot.createContactStatus();

  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.setContactStatusUniformly(contact_status_initial);

  const double t_initial_front_swing = 0.125;
  const double t_initial_front_hip_swing = 0.05;
  const double t_initial_hip_swing = 0.125;
  const double t_initial = t_initial_front_swing + t_initial_front_hip_swing + t_initial_hip_swing;
  const double t_initial_front_swing2 = 0.135;
  const double t_initial_front_hip_swing2 = 0.055;
  const double t_initial_hip_swing2 = 0.15;
  const double t_initial2 = t_initial_front_swing2 + t_initial_front_hip_swing2 + t_initial_hip_swing2;

  contact_status_front_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_swing, t_start, t);
  contact_status_front_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, t_start+t_initial_front_swing, t);

  contact_points[0].coeffRef(0) += 0.25 * stride;
  contact_points[1].coeffRef(0) += 0.25 * stride + 0.5 * additive_stride_hip;
  contact_points[2].coeffRef(0) += 0.25 * stride;
  contact_points[3].coeffRef(0) += 0.25 * stride + 0.5 * additive_stride_hip;

  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, t_start+t_initial_front_swing+t_initial_front_hip_swing, t);

  contact_status_front_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_swing, t_start+t_initial, t);
  contact_status_front_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, t_start+t_initial+t_initial_front_swing2, t);

  contact_points[0].coeffRef(0) += 0.5 * stride;
  contact_points[1].coeffRef(0) += 0.5 * stride + 0.5 * additive_stride_hip;
  contact_points[2].coeffRef(0) += 0.5 * stride;
  contact_points[3].coeffRef(0) += 0.5 * stride + 0.5 * additive_stride_hip;

  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, t_start+t_initial+t_initial_front_swing2+t_initial_front_hip_swing2, t);
  const double t_end_init = t_start+t_initial+t_initial2;

  for (int i=0; i<steps; ++i) {
    contact_status_front_swing.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_front_swing, 
                                     t_end_init+i*t_period, t);
    ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, 
                                     t_end_init+i*t_period+t_front_swing, t);
    contact_points[0].coeffRef(0) += stride;
    contact_points[2].coeffRef(0) += stride;
    contact_points[1].coeffRef(0) += stride;
    contact_points[3].coeffRef(0) += stride;
    contact_status_hip_swing.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_hip_swing, 
                                     t_end_init+i*t_period+t_front_swing+t_front_hip_swing, t);
  }

  contact_status_front_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_swing, 
                                   t_end_init+steps*t_period, t);

  // For the last step
  const double t_end_front_swing = 0.15;
  const double t_end_front_hip_swing = 0.05;
  const double t_end_hip_swing = 0.15;
  const double t_end = t_end_front_swing + t_end_front_hip_swing + t_end_hip_swing;

  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, 
                                   t_end_init+steps*t_period+t_end_front_swing, t);

  contact_points[0].coeffRef(0) += 0.5 * stride;
  contact_points[2].coeffRef(0) += 0.5 * stride;
  contact_points[1].coeffRef(0) += 0.5 * stride - additive_stride_hip;
  contact_points[3].coeffRef(0) += 0.5 * stride - additive_stride_hip;
  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, 
                                   t_end_init+steps*t_period+t_end_front_swing+t_end_front_hip_swing, t);
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_initial, 
                                   t_end_init+steps*t_period+t_end, t);

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << -3, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, // LF
       -0.1, -0.7,  1.0, // LH
        0.1,  0.7, -1.0, // RF
        0.1, -0.7,  1.0; // RH
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp_solver.setStateTrajectory(t, q, v);
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, 350, false);


#ifdef ENABLE_VIEWER
  if (argc != 2) {
    std::cout << "Invalid argment!" << std::endl;
    std::cout << "Package serach path must be specified as the second argment!" << std::endl;
  }
  else {
    const std::string pkg_search_path = argv[1];
    idocp::TrajectoryViewer viewer(pkg_search_path, path_to_urdf);
    const double dt = T/N;
    viewer.display(ocp_solver.getSolution("q"), dt);
  }
#endif 

  return 0;
}