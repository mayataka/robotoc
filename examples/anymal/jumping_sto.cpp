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
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"

#include "robotoc/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "robotoc/utils/trajectory_viewer.hpp"
#endif 


int main(int argc, char *argv[]) {
  const int LF_foot_id = 12;
  const int LH_foot_id = 22;
  const int RF_foot_id = 32;
  const int RH_foot_id = 42;
  std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const double baumgarte_time_step = 0.05;
  robotoc::Robot robot(path_to_urdf, robotoc::BaseJointType::FloatingBase, 
                       contact_frames, baumgarte_time_step);

  const double dt = 0.02;
  const double jump_length = 0.8;
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
  q_ref.coeffRef(0) += jump_length;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 1.0, 0, 0, 1.0, 1.0, 1.0, 
              0.001, 0.001, 0.001, 
              0.001, 0.001, 0.001,
              0.001, 0.001, 0.001,
              0.001, 0.001, 0.001;
  Eigen::VectorXd v_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0);
  Eigen::VectorXd u_weight = Eigen::VectorXd::Constant(robot.dimu(), 1.0e-06);
  Eigen::VectorXd a_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0e-06);
  Eigen::VectorXd qi_weight(Eigen::VectorXd::Zero(robot.dimv()));
  qi_weight << 0, 0, 0, 100.0, 100.0, 100.0,  
               0.1, 0.1, 0.1, 
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1;
  Eigen::VectorXd vi_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0);
  Eigen::VectorXd dvi_weight = Eigen::VectorXd::Constant(robot.dimv(), 1.0e-06);
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_ref(q_ref);
  config_cost->set_q_weight(q_weight);
  config_cost->set_qf_weight(q_weight);
  config_cost->set_qi_weight(qi_weight);
  config_cost->set_v_weight(v_weight);
  config_cost->set_vf_weight(v_weight);
  config_cost->set_vi_weight(vi_weight);
  config_cost->set_dvi_weight(dvi_weight);
  config_cost->set_a_weight(a_weight);
  cost->push_back(config_cost);


  // Create the constraints
  auto constraints           = std::make_shared<robotoc::Constraints>();
  auto joint_position_lower  = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper  = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower  = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper  = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower   = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper   = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone         = std::make_shared<robotoc::FrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);
  constraints->setBarrier(1.0e-03);


  // Create the contact sequence
  const int max_num_impulses = 1;
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot, max_num_impulses);

  robot.updateFrameKinematics(q_standing);
  const Eigen::Vector3d q0_3d_LF = robot.framePosition(LF_foot_id);
  const Eigen::Vector3d q0_3d_LH = robot.framePosition(LH_foot_id);
  const Eigen::Vector3d q0_3d_RF = robot.framePosition(RF_foot_id);
  const Eigen::Vector3d q0_3d_RH = robot.framePosition(RH_foot_id);

  std::vector<Eigen::Vector3d> contact_points = {q0_3d_LF, q0_3d_LH, q0_3d_RF, q0_3d_RH};
  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts({0, 1, 2, 3});
  contact_status_standing.setContactPoints(contact_points);
  contact_sequence->initContactSequence(contact_status_standing);

  auto contact_status_flying = robot.createContactStatus();
  contact_sequence->push_back(contact_status_flying, t0+ground_time-0.3, true);

  contact_points[0].coeffRef(0) += jump_length;
  contact_points[1].coeffRef(0) += jump_length;
  contact_points[2].coeffRef(0) += jump_length;
  contact_points[3].coeffRef(0) += jump_length;
  contact_status_standing.setContactPoints(contact_points);
  contact_sequence->push_back(contact_status_standing, 
                              t0+ground_time+flying_time-0.1, true);

  // Create the STO cost function
  auto sto_cost = std::make_shared<robotoc::STOCostFunction>();
  // Create the STO constraints 
  auto sto_constraints = std::make_shared<robotoc::STOConstraints>(2*max_num_impulses);
  sto_constraints->setMinimumDwellTimes({0.15, 0.15, 0.65});
  sto_constraints->setBarrier(1.0e-03);

  // you can check the contact sequence via
  // std::cout << contact_sequence << std::endl;

  const double T = t0 + flying_time + 2 * ground_time; 
  const int N = std::floor(T / dt);
  const int nthreads = 4;
  robotoc::OCP ocp(robot, contact_sequence, cost, constraints, 
                   sto_cost, sto_constraints, T, N);
  robotoc::OCPSolver ocp_solver(ocp, contact_sequence, nthreads);

  const double t = 0;
  Eigen::VectorXd q(q_standing);
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  ocp_solver.meshRefinement(t);
  ocp_solver.initConstraints(t);

  const bool line_search = false;
  const int num_iteration = 150;
  const double dt_tol_mesh = 0.02;
  const double kkt_tol_mesh = 0.1;
  robotoc::benchmark::convergence(ocp_solver, t, q, v, num_iteration, 
                                  dt_tol_mesh, kkt_tol_mesh, line_search);

  // robotoc::benchmark::CPUTime(ocp_solver, t, q, v, 10000, false);

#ifdef ENABLE_VIEWER
  robotoc::TrajectoryViewer viewer(path_to_urdf, robotoc::BaseJointType::FloatingBase);
  const auto ocp_discretization = ocp_solver.getOCPDiscretization();
  const auto time_steps = ocp_discretization.timeSteps();
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), time_steps, mu);
#endif 

  return 0;
}