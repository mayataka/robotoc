#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/solver/ocp_solver.hpp"
#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/time_varying_task_space_3d_cost.hpp"
#include "idocp/cost/time_varying_com_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"

#include "idocp/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "idocp/utils/trajectory_viewer.hpp"
#endif 

#include "periodic_foot_track_ref.hpp"
#include "periodic_com_ref.hpp"


int main(int argc, char *argv[]) {
  const int LF_foot_id = 12;
  const int LH_foot_id = 22;
  const int RF_foot_id = 32;
  const int RH_foot_id = 42;
  std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const double baumgarte_time_step = 0.04;
  idocp::Robot robot(path_to_urdf, idocp::BaseJointType::FloatingBase, 
                     contact_frames, baumgarte_time_step);

  const double dt = 0.02;
  const double step_length = 0.275;
  const double step_height = 0.125;
  const double period_swing = 0.26;
  const double period_double_support = 0.04;
  const double t0 = 0.10;
  const int cycle = 3; 

  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << 0, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 0, 0, 0, 250000, 250000, 250000, 
              0.0001, 0.0001, 0.0001, 
              0.0001, 0.0001, 0.0001,
              0.0001, 0.0001, 0.0001,
              0.0001, 0.0001, 0.0001;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 100, 100, 100, 100, 100, 100, 
              1, 1, 1, 
              1, 1, 1,
              1, 1, 1,
              1, 1, 1;
  Eigen::VectorXd u_weight = Eigen::VectorXd::Constant(robot.dimu(), 1e-01);
  Eigen::VectorXd qi_weight(Eigen::VectorXd::Zero(robot.dimv()));
  qi_weight << 1, 1, 1, 1, 1, 1,  
               100, 100, 100, 
               100, 100, 100,
               100, 100, 100,
               100, 100, 100;
  Eigen::VectorXd vi_weight = Eigen::VectorXd::Constant(robot.dimv(), 100);
  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  config_cost->set_q_ref(q_standing);
  config_cost->set_q_weight(q_weight);
  config_cost->set_qf_weight(q_weight);
  config_cost->set_qi_weight(qi_weight);
  config_cost->set_v_weight(v_weight);
  config_cost->set_vf_weight(v_weight);
  config_cost->set_vi_weight(vi_weight);
  config_cost->set_u_weight(u_weight);
  cost->push_back(config_cost);

  robot.updateFrameKinematics(q_standing);
  const Eigen::Vector3d q0_3d_LF = robot.framePosition(LF_foot_id);
  const Eigen::Vector3d q0_3d_LH = robot.framePosition(LH_foot_id);
  const Eigen::Vector3d q0_3d_RF = robot.framePosition(RF_foot_id);
  const Eigen::Vector3d q0_3d_RH = robot.framePosition(RH_foot_id);
  const double LF_t0 = t0 + period_swing + period_double_support;
  const double LH_t0 = t0;
  const double RF_t0 = t0 + period_swing + period_double_support;
  const double RH_t0 = t0;
  auto LF_foot_ref = std::make_shared<idocp::PeriodicFootTrackRef>(q0_3d_LF, step_length, step_height, 
                                                                   LF_t0, period_swing, 
                                                                   period_swing+2*period_double_support, false);
  auto LH_foot_ref = std::make_shared<idocp::PeriodicFootTrackRef>(q0_3d_LH, step_length, step_height, 
                                                                   LH_t0, period_swing, 
                                                                   period_swing+2*period_double_support, false);
  auto RF_foot_ref = std::make_shared<idocp::PeriodicFootTrackRef>(q0_3d_RF, step_length, step_height, 
                                                                   RF_t0, period_swing, 
                                                                   period_swing+2*period_double_support, false);
  auto RH_foot_ref = std::make_shared<idocp::PeriodicFootTrackRef>(q0_3d_RH, step_length, step_height, 
                                                                   RH_t0, period_swing, 
                                                                   period_swing+2*period_double_support, false);
  auto LF_cost = std::make_shared<idocp::TimeVaryingTaskSpace3DCost>(robot, LF_foot_id, LF_foot_ref);
  auto LH_cost = std::make_shared<idocp::TimeVaryingTaskSpace3DCost>(robot, LH_foot_id, LH_foot_ref);
  auto RF_cost = std::make_shared<idocp::TimeVaryingTaskSpace3DCost>(robot, RF_foot_id, RF_foot_ref);
  auto RH_cost = std::make_shared<idocp::TimeVaryingTaskSpace3DCost>(robot, RH_foot_id, RH_foot_ref);
  const Eigen::Vector3d foot_track_weight = Eigen::Vector3d::Constant(1.0e06);
  LF_cost->set_q_weight(foot_track_weight);
  LH_cost->set_q_weight(foot_track_weight);
  RF_cost->set_q_weight(foot_track_weight);
  RH_cost->set_q_weight(foot_track_weight);
  cost->push_back(LF_cost);
  cost->push_back(LH_cost);
  cost->push_back(RF_cost);
  cost->push_back(RH_cost);

  Eigen::Vector3d CoM_ref0 = (q0_3d_LF + q0_3d_LH + q0_3d_RF + q0_3d_RH) / 4;
  CoM_ref0(2) = robot.CoM()(2);
  Eigen::Vector3d v_CoM_ref = Eigen::Vector3d::Zero();
  v_CoM_ref.coeffRef(0) = 0.5 * step_length / period_swing;
  auto com_ref = std::make_shared<idocp::PeriodicCoMRef>(CoM_ref0, v_CoM_ref, t0, period_swing, 
                                                         period_double_support, false);
  auto com_cost = std::make_shared<idocp::TimeVaryingCoMCost>(robot, com_ref);
  com_cost->set_q_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost);

  auto constraints           = std::make_shared<idocp::Constraints>();
  auto joint_position_lower  = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper  = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower  = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper  = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower   = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper   = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone         = std::make_shared<idocp::FrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);
  constraints->setBarrier(1.0e-01);

  const double T = t0 + cycle*(2*period_double_support+2*period_swing);
  const int N = T / dt; 
  const int max_num_impulse_phase = 2*cycle;

  const int nthreads = 4;
  const double t = 0;
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, max_num_impulse_phase, nthreads);

  std::vector<Eigen::Vector3d> contact_points = {q0_3d_LF, q0_3d_LH, q0_3d_RF, q0_3d_RH};
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.setContactStatusUniformly(contact_status_initial);

  auto contact_status_even = robot.createContactStatus();
  contact_status_even.activateContacts({0, 2});
  contact_status_even.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_even, t0);

  contact_points[1].coeffRef(0) += step_length;
  contact_points[3].coeffRef(0) += step_length;
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_initial, t0+period_swing);

  auto contact_status_odd = robot.createContactStatus();
  contact_status_odd.activateContacts({1, 3});
  contact_status_odd.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_odd, t0+period_swing+period_double_support);

  contact_points[0].coeffRef(0) += step_length;
  contact_points[2].coeffRef(0) += step_length;
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_initial, t0+2*period_swing+period_double_support);

  for (int i=1; i<cycle; ++i) {
    const double t1 = t0 + i*(2*period_swing+2*period_double_support);
    contact_status_even.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_even, t1);

    contact_points[1].coeffRef(0) += step_length;
    contact_points[3].coeffRef(0) += step_length;
    contact_status_initial.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_initial, t1+period_swing);

    contact_status_odd.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_odd, t1+period_swing+period_double_support);

    contact_points[0].coeffRef(0) += step_length;
    contact_points[2].coeffRef(0) += step_length;
    contact_status_initial.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_initial, t1+2*period_swing+period_double_support);
  }

  Eigen::VectorXd q(q_standing);
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  ocp_solver.initConstraints(t);

  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, 60, line_search);

#ifdef ENABLE_VIEWER
  idocp::TrajectoryViewer viewer(path_to_urdf, idocp::BaseJointType::FloatingBase);
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), dt, mu);
#endif 

  return 0;
}