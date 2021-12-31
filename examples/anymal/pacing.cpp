#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"
#include "robotoc/cost/time_varying_com_cost.hpp"
#include "robotoc/cost/periodic_foot_track_ref.hpp"
#include "robotoc/cost/periodic_com_ref.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/solver/solver_options.hpp"

#include "robotoc/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "robotoc/utils/trajectory_viewer.hpp"
#endif 


int main(int argc, char *argv[]) {
  const int LF_foot_id = 12;
  const int LH_foot_id = 22;
  const int RF_foot_id = 32;
  const int RH_foot_id = 42;
  const std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; 
  const std::vector<robotoc::ContactType> contact_types = {robotoc::ContactType::PointContact, 
                                                           robotoc::ContactType::PointContact,
                                                           robotoc::ContactType::PointContact,
                                                           robotoc::ContactType::PointContact};
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const double baumgarte_time_step = 0.04;
  robotoc::Robot robot(path_to_urdf, robotoc::BaseJointType::FloatingBase, 
                       contact_frames, contact_types, baumgarte_time_step);

  const double dt = 0.02;
  const double step_length = 0.25;
  const double step_height = 0.1;
  const double swing_time = 0.24;
  const double double_support_time = 0.04;
  const double t0 = 0.10;
  const int cycle = 3; 

  // Create the cost function
  auto cost = std::make_shared<robotoc::CostFunction>();
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
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
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
  const Eigen::Vector3d x3d0_LF = robot.framePosition(LF_foot_id);
  const Eigen::Vector3d x3d0_LH = robot.framePosition(LH_foot_id);
  const Eigen::Vector3d x3d0_RF = robot.framePosition(RF_foot_id);
  const Eigen::Vector3d x3d0_RH = robot.framePosition(RH_foot_id);
  const double LF_t0 = t0 + swing_time + double_support_time;
  const double LH_t0 = t0 + swing_time + double_support_time;
  const double RF_t0 = t0;
  const double RH_t0 = t0;
  auto LF_foot_ref = std::make_shared<robotoc::PeriodicFootTrackRef>(x3d0_LF, step_length, step_height, 
                                                                     LF_t0, swing_time, 
                                                                     swing_time+2*double_support_time, false);
  auto LH_foot_ref = std::make_shared<robotoc::PeriodicFootTrackRef>(x3d0_LH, step_length, step_height, 
                                                                     LH_t0, swing_time, 
                                                                     swing_time+2*double_support_time, false);
  auto RF_foot_ref = std::make_shared<robotoc::PeriodicFootTrackRef>(x3d0_RF, step_length, step_height, 
                                                                     RF_t0, swing_time, 
                                                                     swing_time+2*double_support_time, true);
  auto RH_foot_ref = std::make_shared<robotoc::PeriodicFootTrackRef>(x3d0_RH, step_length, step_height, 
                                                                     RH_t0, swing_time, 
                                                                     swing_time+2*double_support_time, true);
  auto LF_cost = std::make_shared<robotoc::TimeVaryingTaskSpace3DCost>(robot, LF_foot_id, LF_foot_ref);
  auto LH_cost = std::make_shared<robotoc::TimeVaryingTaskSpace3DCost>(robot, LH_foot_id, LH_foot_ref);
  auto RF_cost = std::make_shared<robotoc::TimeVaryingTaskSpace3DCost>(robot, RF_foot_id, RF_foot_ref);
  auto RH_cost = std::make_shared<robotoc::TimeVaryingTaskSpace3DCost>(robot, RH_foot_id, RH_foot_ref);
  const Eigen::Vector3d foot_track_weight = Eigen::Vector3d::Constant(1.0e06);
  LF_cost->set_x3d_weight(foot_track_weight);
  LH_cost->set_x3d_weight(foot_track_weight);
  RF_cost->set_x3d_weight(foot_track_weight);
  RH_cost->set_x3d_weight(foot_track_weight);
  cost->push_back(LF_cost);
  cost->push_back(LH_cost);
  cost->push_back(RF_cost);
  cost->push_back(RH_cost);

  Eigen::Vector3d com_ref0 = (x3d0_LF + x3d0_LH + x3d0_RF + x3d0_RH) / 4;
  com_ref0(2) = robot.CoM()(2);
  Eigen::Vector3d vcom_ref = Eigen::Vector3d::Zero();
  vcom_ref.coeffRef(0) = 0.5 * step_length / swing_time;
  auto com_ref = std::make_shared<robotoc::PeriodicCoMRef>(com_ref0, vcom_ref, t0, swing_time, 
                                                           double_support_time, true);
  auto com_cost = std::make_shared<robotoc::TimeVaryingCoMCost>(robot, com_ref);
  com_cost->set_com_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost);

  // Create the constraints
  const double barrier = 1.0e-03;
  const double fraction_to_boundary_rule = 0.995;
  auto constraints          = std::make_shared<robotoc::Constraints>(barrier, fraction_to_boundary_rule);
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone        = std::make_shared<robotoc::FrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);

  // Create the contact sequence
  const int max_num_impulses = 2*cycle;
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot, max_num_impulses);

  std::vector<Eigen::Vector3d> contact_positions = {x3d0_LF, x3d0_LH, x3d0_RF, x3d0_RH};
  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts({0, 1, 2, 3});
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->initContactSequence(contact_status_standing);

  auto contact_status_rfrh_swing = robot.createContactStatus();
  contact_status_rfrh_swing.activateContacts({0, 1});
  contact_status_rfrh_swing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_rfrh_swing, t0);

  contact_positions[2].coeffRef(0) += 0.5 * step_length;
  contact_positions[3].coeffRef(0) += 0.5 * step_length;
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_standing, t0+swing_time);

  auto contact_status_lflh_swing = robot.createContactStatus();
  contact_status_lflh_swing.activateContacts({2, 3});
  contact_status_lflh_swing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_lflh_swing, 
                              t0+swing_time+double_support_time);

  contact_positions[0].coeffRef(0) += step_length;
  contact_positions[1].coeffRef(0) += step_length;
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_standing, 
                              t0+2*swing_time+double_support_time);

  for (int i=1; i<cycle; ++i) {
    const double t1 = t0 + i*(2*swing_time+2*double_support_time);
    contact_status_rfrh_swing.setContactPlacements(contact_positions);
    contact_sequence->push_back(contact_status_rfrh_swing, t1);

    contact_positions[2].coeffRef(0) += step_length;
    contact_positions[3].coeffRef(0) += step_length;
    contact_status_standing.setContactPlacements(contact_positions);
    contact_sequence->push_back(contact_status_standing, t1+swing_time);

    contact_status_lflh_swing.setContactPlacements(contact_positions);
    contact_sequence->push_back(contact_status_lflh_swing, 
                                t1+swing_time+double_support_time);

    contact_positions[0].coeffRef(0) += step_length;
    contact_positions[1].coeffRef(0) += step_length;
    contact_status_standing.setContactPlacements(contact_positions);
    contact_sequence->push_back(contact_status_standing, 
                                t1+2*swing_time+double_support_time);
  }

  // you can check the contact sequence via
  // std::cout << contact_sequence << std::endl;

  // Create the OCP solver.
  const double T = t0 + cycle*(2*double_support_time+2*swing_time);
  const int N = T / dt; 
  robotoc::OCP ocp(robot, cost, constraints, T, N, max_num_impulses);
  auto solver_options = robotoc::SolverOptions::defaultOptions();
  const int nthreads = 4;
  robotoc::OCPSolver ocp_solver(ocp, contact_sequence, solver_options, nthreads);

  // Initial time and initial state
  const double t = 0;
  const Eigen::VectorXd q(q_standing);
  const Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  // Solves the OCP.
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);
  ocp_solver.initConstraints(t);
  std::cout << "Initial KKT error: " << ocp_solver.KKTError(t, q, v) << std::endl;
  ocp_solver.solve(t, q, v);
  std::cout << "KKT error after convergence: " << ocp_solver.KKTError(t, q, v) << std::endl;
  std::cout << ocp_solver.getSolverStatistics() << std::endl;

  // const int num_iteration = 10000;
  // robotoc::benchmark::CPUTime(ocp_solver, t, q, v, num_iteration);

#ifdef ENABLE_VIEWER
  robotoc::TrajectoryViewer viewer(path_to_urdf, robotoc::BaseJointType::FloatingBase);
  const auto discretization = ocp_solver.getTimeDiscretization();
  const auto time_steps = discretization.timeSteps();
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), time_steps, mu);
#endif 

  return 0;
}