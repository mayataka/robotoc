#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
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
  const int LF_foot_id = 14;
  const int LH_foot_id = 24;
  const int RF_foot_id = 34;
  const int RH_foot_id = 44;
  std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const double baumgarte_time_step = 0.04;
  idocp::Robot robot(path_to_urdf, contact_frames, baumgarte_time_step);

  const double dt = 0.01;
  const double jump_length = 0.5;
  const double jump_height = 0.1;
  const double period_flying_up = 0.15;
  const double period_flying_down = period_flying_up;
  const double period_flying = period_flying_up + period_flying_down;
  const double period_ground = 0.30;
  const double t0 = 0;

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

  Eigen::Vector3d CoM_ref0_flying_up = (q0_3d_LF + q0_3d_LH + q0_3d_RF + q0_3d_RH) / 4;
  CoM_ref0_flying_up(2) = robot.CoM()(2);
  Eigen::Vector3d v_CoM_ref_flying_up = Eigen::Vector3d::Zero();
  v_CoM_ref_flying_up << (0.5*jump_length/period_flying_up), 0, (jump_height/period_flying_up);
  auto com_ref_flying_up = std::make_shared<idocp::PeriodicCoMRef>(CoM_ref0_flying_up, v_CoM_ref_flying_up, 
                                                                   t0+period_ground, period_flying_up, 
                                                                   period_flying_down+2*period_ground, false);
  auto com_cost_flying_up = std::make_shared<idocp::TimeVaryingCoMCost>(robot, com_ref_flying_up);
  com_cost_flying_up->set_q_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost_flying_up);

  Eigen::Vector3d CoM_ref0_landed = (q0_3d_LF + q0_3d_LH + q0_3d_RF + q0_3d_RH) / 4;
  CoM_ref0_landed(0) += jump_length;
  CoM_ref0_landed(2) = robot.CoM()(2);
  const Eigen::Vector3d v_CoM_ref_landed = Eigen::Vector3d::Zero();
  auto com_ref_landed = std::make_shared<idocp::PeriodicCoMRef>(CoM_ref0_landed, v_CoM_ref_landed, 
                                                                t0+period_ground+period_flying, period_ground, 
                                                                period_ground+period_flying, false);
  auto com_cost_landed = std::make_shared<idocp::TimeVaryingCoMCost>(robot, com_ref_landed);
  com_cost_landed->set_q_weight(Eigen::Vector3d::Constant(1.0e06));
  cost->push_back(com_cost_landed);

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
  constraints->setBarrier(1.0e-03);

  const double T = t0 + period_flying + 2 * period_ground; 
  const int N = T / dt;
  const int max_num_impulse_phase = 1;

  const int nthreads = 4;
  const double t = 0;
  idocp::OCPSolver ocp_solver(robot, cost, constraints, T, N, max_num_impulse_phase, nthreads);

  std::vector<Eigen::Vector3d> contact_points = {q0_3d_LF, q0_3d_LH, q0_3d_RF, q0_3d_RH};
  auto contact_status_initial = robot.createContactStatus();
  contact_status_initial.activateContacts({0, 1, 2, 3});
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.setContactStatusUniformly(contact_status_initial);

  auto contact_status_flying = robot.createContactStatus();
  ocp_solver.pushBackContactStatus(contact_status_flying, t0+period_ground);

  contact_points[0].coeffRef(0) += jump_length;
  contact_points[1].coeffRef(0) += jump_length;
  contact_points[2].coeffRef(0) += jump_length;
  contact_points[3].coeffRef(0) += jump_length;
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_initial, t0+period_ground+period_flying);

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
  idocp::TrajectoryViewer viewer(path_to_urdf);
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), dt, mu);
#endif 

  return 0;
}