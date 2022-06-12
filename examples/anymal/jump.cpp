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
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const std::vector<std::string> contact_frames = {"LF_FOOT", "LH_FOOT", "RF_FOOT", "RH_FOOT"}; 
  const std::vector<robotoc::ContactType> contact_types = {robotoc::ContactType::PointContact, 
                                                           robotoc::ContactType::PointContact,
                                                           robotoc::ContactType::PointContact,
                                                           robotoc::ContactType::PointContact};
  const double baumgarte_time_step = 0.04;
  robotoc::Robot robot(path_to_urdf, robotoc::BaseJointType::FloatingBase, 
                       contact_frames, contact_types, baumgarte_time_step);

  const double dt = 0.01;
  const Eigen::Vector3d jump_length = {0.5, 0, 0};
  const double jump_height = 0.1;
  const double flying_up_time = 0.15;
  const double flying_down_time = flying_up_time;
  const double flying_time = flying_up_time + flying_down_time;
  const double ground_time = 0.30;
  const double t0 = 0;

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
  const Eigen::Vector3d x3d0_LF = robot.framePosition("LF_FOOT");
  const Eigen::Vector3d x3d0_LH = robot.framePosition("LH_FOOT");
  const Eigen::Vector3d x3d0_RF = robot.framePosition("RF_FOOT");
  const Eigen::Vector3d x3d0_RH = robot.framePosition("RH_FOOT");

  const Eigen::Vector3d com_ref0_flying_up = robot.CoM();
  const Eigen::Vector3d vcom_ref_flying_up = 0.5*jump_length/flying_up_time 
                                            + Eigen::Vector3d({0, 0, jump_height/flying_up_time});
  auto com_ref_flying_up = std::make_shared<robotoc::PeriodicCoMRef>(com_ref0_flying_up, vcom_ref_flying_up, 
                                                                     t0+ground_time, flying_up_time, 
                                                                     flying_down_time+2*ground_time, false);
  auto com_cost_flying_up = std::make_shared<robotoc::TimeVaryingCoMCost>(robot, com_ref_flying_up);
  com_cost_flying_up->set_com_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost_flying_up);

  const Eigen::Vector3d com_ref0_landed = robot.CoM() + jump_length;
  const Eigen::Vector3d vcom_ref_landed = Eigen::Vector3d::Zero();
  auto com_ref_landed = std::make_shared<robotoc::PeriodicCoMRef>(com_ref0_landed, vcom_ref_landed, 
                                                                  t0+ground_time+flying_time, ground_time, 
                                                                  ground_time+flying_time, false);
  auto com_cost_landed = std::make_shared<robotoc::TimeVaryingCoMCost>(robot, com_ref_landed);
  com_cost_landed->set_com_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost_landed);

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
  const int max_num_each_discrete_events = 1;
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot, max_num_each_discrete_events);

  std::unordered_map<std::string, Eigen::Vector3d> contact_positions 
      = {{"LF_FOOT", x3d0_LF}, {"LH_FOOT", x3d0_LH}, {"RF_FOOT", x3d0_RF}, {"RH_FOOT", x3d0_RH}};
  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts(std::vector<std::string>({"LF_FOOT", "LH_FOOT", "RF_FOOT", "RH_FOOT"}));
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->initContactSequence(contact_status_standing);

  auto contact_status_flying = robot.createContactStatus();
  contact_sequence->push_back(contact_status_flying, t0+ground_time);

  contact_positions["LF_FOOT"].noalias() += jump_length;
  contact_positions["LH_FOOT"].noalias() += jump_length;
  contact_positions["RF_FOOT"].noalias() += jump_length;
  contact_positions["RH_FOOT"].noalias() += jump_length;
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_standing, 
                              t0+ground_time+flying_time);

  // you can check the contact sequence via
  // std::cout << contact_sequence << std::endl;

  // Create the OCP solver.
  const double T = t0 + flying_time + 2 * ground_time; 
  const int N = T / dt;
  robotoc::OCP ocp(robot, cost, constraints, T, N, max_num_each_discrete_events);
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