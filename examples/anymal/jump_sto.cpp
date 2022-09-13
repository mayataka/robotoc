#include <string>
#include <memory>

#include "Eigen/Core"

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/task_space_3d_cost.hpp"
#include "robotoc/cost/com_cost.hpp"
#include "robotoc/cost/periodic_swing_foot_ref.hpp"
#include "robotoc/cost/periodic_com_ref.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"
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
  const double baumgarte_time_step = 0.05;
  robotoc::Robot robot(path_to_urdf, robotoc::BaseJointType::FloatingBase, 
                       contact_frames, contact_types, baumgarte_time_step);

  const double dt = 0.02;
  const Eigen::Vector3d jump_length = {0.8, 0, 0};
  const double flying_up_time = 0.15;
  const double flying_down_time = flying_up_time;
  const double flying_time = flying_up_time + flying_down_time;
  const double ground_time = 0.70;
  const double t0 = 0;

  // Create the cost function
  auto cost = std::make_shared<robotoc::CostFunction>();
  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << 0, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_ref = q_standing;
  q_ref.head(3).noalias() += jump_length;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 1.0, 0, 0, 1.0, 1.0, 1.0, 
              0.001, 0.001, 0.001, 
              0.001, 0.001, 0.001,
              0.001, 0.001, 0.001,
              0.001, 0.001, 0.001;
  Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0);
  Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0e-06);
  Eigen::VectorXd q_weight_impulse(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight_impulse << 0, 0, 0, 100.0, 100.0, 100.0,  
               0.1, 0.1, 0.1, 
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1;
  Eigen::VectorXd v_weight_impulse = Eigen::VectorXd::Constant(robot.dimv(), 1.0);
  Eigen::VectorXd dv_weight_impulse = Eigen::VectorXd::Constant(robot.dimv(), 1.0e-06);
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_ref(q_ref);
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_weight_terminal(q_weight);
  config_cost->set_q_weight_impulse(q_weight_impulse);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_weight_terminal(v_weight);
  config_cost->set_v_weight_impulse(v_weight_impulse);
  config_cost->set_dv_weight_impulse(dv_weight_impulse);
  config_cost->set_a_weight(a_weight);
  cost->push_back(config_cost);

  // Create the constraints
  const double barrier_param = 1.0e-03;
  const double fraction_to_boundary_rule = 0.995;
  auto constraints          = std::make_shared<robotoc::Constraints>(barrier_param, fraction_to_boundary_rule);
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  auto friction_cone        = std::make_shared<robotoc::FrictionCone>(robot);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);

  // Create the contact sequence
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot);
  const double mu = 0.7;
  const std::unordered_map<std::string, double> friction_coefficients = {{"LF_FOOT", mu}, 
                                                                         {"LH_FOOT", mu}, 
                                                                         {"RF_FOOT", mu}, 
                                                                         {"RH_FOOT", mu}};

  robot.updateFrameKinematics(q_standing);
  const Eigen::Vector3d x3d0_LF = robot.framePosition("LF_FOOT");
  const Eigen::Vector3d x3d0_LH = robot.framePosition("LH_FOOT");
  const Eigen::Vector3d x3d0_RF = robot.framePosition("RF_FOOT");
  const Eigen::Vector3d x3d0_RH = robot.framePosition("RH_FOOT");

  std::unordered_map<std::string, Eigen::Vector3d> contact_positions = {{"LF_FOOT", x3d0_LF}, 
                                                                        {"LH_FOOT", x3d0_LH}, 
                                                                        {"RF_FOOT", x3d0_RF}, 
                                                                        {"RH_FOOT", x3d0_RH}};
  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts(std::vector<std::string>({"LF_FOOT", "LH_FOOT", "RF_FOOT", "RH_FOOT"}));
  contact_status_standing.setContactPlacements(contact_positions);
  contact_status_standing.setFrictionCoefficients(friction_coefficients);
  contact_sequence->init(contact_status_standing);

  auto contact_status_flying = robot.createContactStatus();
  contact_sequence->push_back(contact_status_flying, t0+ground_time-0.3, true);

  contact_positions["LF_FOOT"].noalias() += jump_length;
  contact_positions["LH_FOOT"].noalias() += jump_length;
  contact_positions["RF_FOOT"].noalias() += jump_length;
  contact_positions["RH_FOOT"].noalias() += jump_length;
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->push_back(contact_status_standing, 
                              t0+ground_time+flying_time-0.1, true);

  // Create the STO cost function
  auto sto_cost = std::make_shared<robotoc::STOCostFunction>();
  // Create the STO constraints 
  const std::vector<double> min_dwell_times = {0.15, 0.15, 0.65};
  auto sto_constraints = std::make_shared<robotoc::STOConstraints>(min_dwell_times,
                                                                   barrier_param, 
                                                                   fraction_to_boundary_rule);

  // you can check the contact sequence via
  // std::cout << contact_sequence << std::endl;

  const double T = t0 + flying_time + 2 * ground_time; 
  const int N = std::floor(T / dt);
  robotoc::OCP ocp(robot, cost, constraints, sto_cost, sto_constraints, 
                   contact_sequence, T, N);
  auto solver_options = robotoc::SolverOptions();
  solver_options.max_dt_mesh = T/N;
  solver_options.kkt_tol_mesh = 0.1;
  solver_options.max_iter = 200;
  const int nthreads = 4;
  robotoc::OCPSolver ocp_solver(ocp, solver_options, nthreads);

  // Initial time and initial state
  const double t = 0;
  Eigen::VectorXd q(q_standing);
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  // Solves the OCP.
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);
  ocp_solver.meshRefinement(t);
  ocp_solver.initConstraints(t);
  std::cout << "Initial KKT error: " << ocp_solver.KKTError(t, q, v) << std::endl;
  ocp_solver.solve(t, q, v);
  std::cout << "KKT error after convergence: " << ocp_solver.KKTError(t, q, v) << std::endl;
  std::cout << ocp_solver.getSolverStatistics() << std::endl;

  // const int num_iteration = 10000;
  // robotoc::benchmark::CPUTime(ocp_solver, t, q, v, num_iteration);

#ifdef ENABLE_VIEWER
  robotoc::TrajectoryViewer viewer(path_to_urdf, robotoc::BaseJointType::FloatingBase);
  const auto time_discretization = ocp_solver.getTimeDiscretization();
  const auto time_steps = time_discretization.timeSteps();
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), time_steps, mu);
#endif 

  return 0;
}