#include <string>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp_solver.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/time_varying_configuration_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"
#include "idocp/constraints/impulse_friction_cone.hpp"

#include "idocp/utils/ocp_benchmarker.hpp"

#ifdef ENABLE_VIEWER
#include "idocp/utils/trajectory_viewer.hpp"
#endif 

 
class TimeVaryingConfigurationRef final : public idocp::TimeVaryingConfigurationRefBase {
public:
  TimeVaryingConfigurationRef(const double t0, const double tf, 
                              const Eigen::VectorXd& q0, const double v_ref) 
    : TimeVaryingConfigurationRefBase(),
      t0_(t0),
      tf_(tf),
      q0_(q0),
      qf_(q0),
      v_ref_(v_ref) {
    qf_.coeffRef(0) += (tf-t0) * v_ref;
  }

  ~TimeVaryingConfigurationRef() {}

  void update_q_ref(const idocp::Robot& robot, const double t, 
                    Eigen::VectorXd& q_ref) const override {
    if (t < t0_) {
      q_ref = q0_;
    }
    else if (t < tf_) {
      q_ref = q0_;
      q_ref.coeffRef(0) += (t-t0_) * v_ref_;
    }
    else {
      q_ref = qf_;
    }
  }

  bool isActive(const double t) const override {
    return true;
  }

private:
  Eigen::VectorXd q0_, qf_;
  double t0_, tf_, v_ref_;
};


int main(int argc, char *argv[]) {
  const int LF_foot_id = 14;
  const int LH_foot_id = 24;
  const int RF_foot_id = 34;
  const int RH_foot_id = 44;
  std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; // LF, LH, RF, RH
  const std::string path_to_urdf = "../anymal_b_simple_description/urdf/anymal.urdf";
  const double baumgarte_time_step = 0.04;
  idocp::Robot robot(path_to_urdf, contact_frames, baumgarte_time_step);

  const double stride = 0.45;
  const double additive_stride_hip = 0.2;
  const double t_start = 1.0;

  const double t_front_swing = 0.135;
  const double t_front_hip_swing = 0.05;
  const double t_hip_swing = 0.165;
  const double t_period = t_front_swing + t_front_hip_swing + t_hip_swing;
  const int steps = 10;

  auto cost = std::make_shared<idocp::CostFunction>();
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.001));

  auto config_cost = std::make_shared<idocp::ConfigurationSpaceCost>(robot);
  config_cost->set_v_weight(v_weight);
  config_cost->set_vf_weight(v_weight);
  config_cost->set_vi_weight(v_weight);
  config_cost->set_a_weight(a_weight);
  config_cost->set_dvi_weight(a_weight);
  cost->push_back(config_cost);

  Eigen::VectorXd q_standing(Eigen::VectorXd::Zero(robot.dimq()));
  q_standing << -3, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 100, 100, 100, 100, 100, 100, 
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1;
  const double v_ref = stride / t_period;
  auto config_ref = std::make_shared<TimeVaryingConfigurationRef>(t_start+0.25*t_period, 0.25*t_period+t_start+(0.75+steps+0.75)*t_period, q_standing, v_ref);
  auto time_varying_config_cost = std::make_shared<idocp::TimeVaryingConfigurationSpaceCost>(robot, config_ref);
  time_varying_config_cost->set_q_weight(q_weight);
  time_varying_config_cost->set_qf_weight(q_weight);
  time_varying_config_cost->set_qi_weight(q_weight);
  cost->push_back(time_varying_config_cost);

  auto constraints           = std::make_shared<idocp::Constraints>();
  auto joint_position_lower  = std::make_shared<idocp::JointPositionLowerLimit>(robot);
  auto joint_position_upper  = std::make_shared<idocp::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower  = std::make_shared<idocp::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper  = std::make_shared<idocp::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower   = std::make_shared<idocp::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper   = std::make_shared<idocp::JointTorquesUpperLimit>(robot);
  const double mu = 0.7;
  auto friction_cone         = std::make_shared<idocp::FrictionCone>(robot, mu);
  auto impulse_friction_cone = std::make_shared<idocp::ImpulseFrictionCone>(robot, mu);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(friction_cone);
  constraints->push_back(impulse_friction_cone);
  constraints->setBarrier(1.0e-03);

  const double T = 7; 
  const int N = 240;
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
  ocp_solver.pushBackContactStatus(contact_status_front_swing, t_start);
  contact_status_front_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, t_start+t_initial_front_swing);

  contact_points[0].coeffRef(0) += 0.25 * stride;
  contact_points[1].coeffRef(0) += 0.25 * stride + 0.5 * additive_stride_hip;
  contact_points[2].coeffRef(0) += 0.25 * stride;
  contact_points[3].coeffRef(0) += 0.25 * stride + 0.5 * additive_stride_hip;

  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, t_start+t_initial_front_swing+t_initial_front_hip_swing);

  contact_status_front_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_swing, t_start+t_initial);
  contact_status_front_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, t_start+t_initial+t_initial_front_swing2);

  contact_points[0].coeffRef(0) += 0.5 * stride;
  contact_points[1].coeffRef(0) += 0.5 * stride + 0.5 * additive_stride_hip;
  contact_points[2].coeffRef(0) += 0.5 * stride;
  contact_points[3].coeffRef(0) += 0.5 * stride + 0.5 * additive_stride_hip;

  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, t_start+t_initial+t_initial_front_swing2+t_initial_front_hip_swing2);
  const double t_end_init = t_start+t_initial+t_initial2;

  for (int i=0; i<steps; ++i) {
    contact_status_front_swing.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_front_swing, t_end_init+i*t_period);
    ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, t_end_init+i*t_period+t_front_swing);
    contact_points[0].coeffRef(0) += stride;
    contact_points[2].coeffRef(0) += stride;
    contact_points[1].coeffRef(0) += stride;
    contact_points[3].coeffRef(0) += stride;
    contact_status_hip_swing.setContactPoints(contact_points);
    ocp_solver.pushBackContactStatus(contact_status_hip_swing, 
                                     t_end_init+i*t_period+t_front_swing+t_front_hip_swing);
  }

  contact_status_front_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_front_swing, t_end_init+steps*t_period);

  // For the last step
  const double t_end_front_swing = 0.15;
  const double t_end_front_hip_swing = 0.05;
  const double t_end_hip_swing = 0.15;
  const double t_end = t_end_front_swing + t_end_front_hip_swing + t_end_hip_swing;

  ocp_solver.pushBackContactStatus(contact_status_front_hip_swing, 
                                   t_end_init+steps*t_period+t_end_front_swing);

  contact_points[0].coeffRef(0) += 0.75 * stride;
  contact_points[2].coeffRef(0) += 0.75 * stride;
  contact_points[1].coeffRef(0) += 0.75 * stride - additive_stride_hip;
  contact_points[3].coeffRef(0) += 0.75 * stride - additive_stride_hip;
  contact_status_hip_swing.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_hip_swing, 
                                   t_end_init+steps*t_period+t_end_front_swing+t_end_front_hip_swing);
  contact_status_initial.setContactPoints(contact_points);
  ocp_solver.pushBackContactStatus(contact_status_initial, t_end_init+steps*t_period+t_end);

  Eigen::VectorXd q(q_standing);
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);
  ocp_solver.setSolution("lmd", f_init);

  ocp_solver.initConstraints(t);

  const bool line_search = false;
  idocp::ocpbenchmarker::Convergence(ocp_solver, t, q, v, 200, line_search);
  // idocp::ocpbenchmarker::CPUTime(ocp_solver, t, q, v, 2500, line_search);

#ifdef ENABLE_VIEWER
  idocp::TrajectoryViewer viewer(path_to_urdf);
  Eigen::Vector3d camera_pos;
  Eigen::Vector4d camera_quat;
  camera_pos << 5.10483, -3.98692, 1.59321;
  camera_quat << 0.547037, 0.243328, 0.314829, 0.736495;
  viewer.setCameraTransform(camera_pos, camera_quat);
  viewer.display(robot, ocp_solver.getSolution("q"), 
                 ocp_solver.getSolution("f", "WORLD"), (T/N), mu);
  camera_pos << 0.119269, -7.96283, 1.95978;
  camera_quat << 0.609016, 0.00297497, 0.010914, 0.793077;
  viewer.setCameraTransform(camera_pos, camera_quat);
  viewer.display(ocp_solver.getSolution("q"), (T/N));
#endif 

  return 0;
}