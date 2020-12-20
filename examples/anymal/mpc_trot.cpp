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
#include "idocp/cost/time_varying_configuration_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/friction_cone.hpp"
#include "idocp/constraints/contact_normal_force.hpp"

#include "idocp/utils/quadruped_simulator.hpp"


class MPCCallbackStanding {
public:
  MPCCallbackStanding(const idocp::Robot& robot)
    : robot_(robot),
      contact_points_(robot.maxPointContacts(), Eigen::Vector3d::Zero()) {} 

  template <typename OCPSolverType>
  void init(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            idocp::MPC<OCPSolverType>& mpc) {}

  template <typename OCPSolverType>
  void callback(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v, idocp::MPC<OCPSolverType>& mpc) {
    robot_.updateFrameKinematics(q);
    robot_.getContactPoints(contact_points_);
    mpc.getSolverHandle()->setContactPoints(0, contact_points_);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  idocp::Robot robot_;
  std::vector<Eigen::Vector3d> contact_points_;
};


class MPCCallbackTrotting {
public:
  MPCCallbackTrotting(const idocp::Robot& robot)
    : robot_(robot),
      contact_points_init_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
      contact_points_even_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
      contact_points_odd_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
      contact_status_init_(robot.createContactStatus()),
      contact_status_even_(robot.createContactStatus()),
      contact_status_odd_(robot.createContactStatus()),
      steps_(0),
      step_time_(),
      contact_frames_({14, 24, 34, 44}) {
    contact_status_even_.activateContacts({0, 1});
    contact_status_odd_.activateContacts({2, 3});
  } 

  template <typename OCPSolverType>
  void init(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            idocp::MPC<OCPSolverType>& mpc) {
    step_time_.push_back(t+t_start);
  }

  template <typename OCPSolverType>
  void callback(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v, idocp::MPC<OCPSolverType>& mpc) {
    robot_.updateFrameKinematics(q);
    robot_.getContactPoints(contact_points_even_);
    mpc.getSolverHandle()->setContactPoints(0, contact_points_even_);
    // if (t <= 0.5) {
    //   mpc.getSolverHandle()->setContactPoints(1, contact_points_even_);
    // }
    if (step_time_.front() <= t) {
      step_time_.pop_front();
      mpc.getSolverHandle()->popFrontDiscreteEvent();
    }
    if (step_time_.back()+kPeriod < t+kT) {
      if (steps_%2 == 0) {
      }
      else {
      }
      // mpc.getSolverHandle()->pushBackDiscreteEvent();
    }
    // if (steps_%2 == 0) {
    //   robot_.updateFrameKinematics(q);
    //   robot_.getContactPoints(contact_points_even_);
    //   contact_status_even_.getContactPoints(contact_points_even_);
    //   mpc.getSolverHandle()->setContactPoints(0, contact_points_even_);
    //   mpc.getSolverHandle()->setContactPoints(1, contact_points_odd_);
    // }
    // else {
    //   robot_.updateFrameKinematics(q);= {14, 24, 34, 44}
    //   robot_.getContactPoints(contact_points_odd_);
    //   mpc.getSolverHandle()->setContactPoints(0, contact_points_odd_);
    //   mpc.getSolverHandle()->setContactPoints(1, contact_points_even_);
    // }
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  idocp::Robot robot_;
  std::vector<Eigen::Vector3d> contact_points_init_, contact_points_even_, contact_points_odd_;
  idocp::ContactStatus contact_status_init_, contact_status_even_, contact_status_odd_;
  int steps_;
  std::deque<double> step_time_;
  static constexpr double kPeriod = 0.5;
  static constexpr double kStepLength = 0.15;
  static constexpr double kT = 0.51;
  static constexpr double t_start = 0.5;
  std::vector<int> contact_frames_;
};



int main(int argc, char *argv[]) {
  srand((unsigned int) time(0));
  std::vector<int> contact_frames = {14, 24, 34, 44};
  const std::string path_to_urdf = "../anymal/anymal.urdf";
  const std::string path_to_urdf_for_raisim = "/home/sotaro/src/idocp/examples/anymal/anymal/anymal_for_raisim.urdf";
  idocp::Robot robot(path_to_urdf, contact_frames);
  auto cost = std::make_shared<idocp::CostFunction>();
  // Eigen::VectorXd q_ref(Eigen::VectorXd::Zero(robot.dimq()));
  // q_ref << 0, 0, 0.4792, 0, 0, 0, 1, 
  //          -0.1,  0.7, -1.0, 
  //          -0.1, -0.7,  1.0, 
  //           0.1,  0.7, -1.0, 
  //           0.1, -0.7,  1.0;
  // Eigen::VectorXd v_ref(Eigen::VectorXd::Zero(robot.dimv()));
  // v_ref << 0.5, 0, 0, 0, 0, 0, 
  //            0,  2.0, -2.25,
  //            0,  1.5,  0.5,
  //            0, -1.5, -0.5,
  //            0, -2.0,  2.25;
  // auto configuration_cost = std::make_shared<idocp::TimeVaryingConfigurationCost>(robot);
  // configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // configuration_cost->set_ref(0, q_ref, v_ref);
  // configuration_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // configuration_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // configuration_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  // cost->push_back(configuration_cost);

  Eigen::VectorXd q_ref(Eigen::VectorXd::Zero(robot.dimq()));
  q_ref << 0, 0, 0.4792, 0, 0, 0, 1, 
           -0.1,  0.7, -1.0, 
           -0.1, -0.7,  1.0, 
            0.1,  0.7, -1.0, 
            0.1, -0.7,  1.0;
  // q_ref << 0.0, 0, 0.4792, 0, 0, 0, 1, 
  //          -0.1,  0.0, -1.0, 
  //          -0.1, -0.7,  1.0, 
  //           0.1,  0.7, -1.0, 
  //           0.1, -1.5,  1.5;
  Eigen::VectorXd v_ref(Eigen::VectorXd::Zero(robot.dimv()));
  v_ref << 0.0, 0, 0, 0, 0, 0, 
             0,  0,  0,
             0,  0,  0,
             0,  0,  0,
             0,  0,  0;
  auto joint_cost = std::make_shared<idocp::JointSpaceCost>(robot);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_ref(v_ref);
  // Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  // q_weight << 10, 10, 10, 10, 10, 10, 
  //            0,  0,  0,
  //            0,  0,  0,
  //            0,  0,  0,
  //            0,  0,  0;
  // joint_cost->set_q_weight(q_weight);
  // joint_cost->set_qf_weight(q_weight);
  joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_qf_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_vf_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  joint_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  joint_cost->set_u_weight(Eigen::VectorXd::Constant(robot.dimu(), 0.0));
  cost->push_back(joint_cost);

  auto contact_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.01, 0.01, 0.01;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    // fr << 0, 0, 0;
    fr << 0, 0, 70;
    f_ref.push_back(fr);
  }
  contact_cost->set_f_weight(f_weight);
  contact_cost->set_f_ref(f_ref);
  cost->push_back(contact_cost);
  auto impulse_joint_cost = std::make_shared<idocp::JointSpaceImpulseCost>(robot);
  impulse_joint_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  impulse_joint_cost->set_q_ref(q_ref);
  impulse_joint_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  impulse_joint_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  cost->push_back(impulse_joint_cost);
  // auto impulse_configuration_cost = std::make_shared<idocp::ImpulseTimeVaryingConfigurationCost>(robot);
  // impulse_configuration_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  // impulse_configuration_cost->set_ref(0, q_ref, v_ref);
  // impulse_configuration_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  // impulse_configuration_cost->set_dv_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.1));
  // cost->push_back(impulse_configuration_cost);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  impulse_force_cost->set_f_ref(f_ref);
  cost->push_back(impulse_force_cost);

  auto task_3d_cost1 = std::make_shared<idocp::TaskSpace3DCost>(robot, 14);
  task_3d_cost1->set_q_3d_weight(Eigen::Vector3d::Constant(100));
  task_3d_cost1->set_qf_3d_weight(Eigen::Vector3d::Constant(100));
  Eigen::Vector3d q_ref1;
  q_ref1 << 0.369915+0.1, 0.198573, 0.01;
  task_3d_cost1->set_q_3d_ref(q_ref1);
  cost->push_back(task_3d_cost1);

  // auto task_3d_cost4 = std::make_shared<idocp::TaskSpace3DCost>(robot, 44);
  // task_3d_cost4->set_q_3d_weight(Eigen::Vector3d::Constant(100));
  // task_3d_cost4->set_qf_3d_weight(Eigen::Vector3d::Constant(100));
  // Eigen::Vector3d q_ref4;
  // q_ref4 << -0.369915+0.15, -0.198573, 0.05;
  // task_3d_cost4->set_q_3d_ref(q_ref4);
  // cost->push_back(task_3d_cost4);

  auto constraints = std::make_shared<idocp::Constraints>();
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
  const double T = 0.6;
  const int N = 20;
  const int max_num_impulse_phase = 4;
  const int num_proc = 4;
  const double t = 0;
  idocp::MPC<idocp::OCPSolver> mpc(robot, cost, constraints, T, N, max_num_impulse_phase, num_proc);
  robot.printRobotModel();

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q <<  0, 0, 0.4792, 0, 0, 0, 1, 
       -0.1,  0.7, -1.0, 
       -0.1, -0.7,  1.0, 
        0.1,  0.7, -1.0, 
        0.1, -0.7,  1.0;
  Eigen::VectorXd v(Eigen::VectorXd::Zero(robot.dimv()));

  auto contact_status = robot.createContactStatus();
  contact_status.activateContacts({0, 1, 2, 3});
  robot.updateFrameKinematics(q);
  std::vector<Eigen::Vector3d> contact_points(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  robot.getContactPoints(contact_points);
  contact_status.setContactPoints(contact_points);
  mpc.getSolverHandle()->setContactStatusUniformly(contact_status);

  // auto contact_status_next = robot.createContactStatus();
  // contact_status_next.setContactPoints(contact_points);
  // contact_status_next.activateContacts({1, 2, 3});
  // constexpr double switchig_time = 0.5;
  // mpc.getSolverHandle()->pushBackContactStatus(contact_status_next, switchig_time, 0);

  mpc.initializeSolution(t, q, v, 100);
  mpc.getSolverHandle()->printSolution("end-effector", {14, 24, 34, 44});
  const std::string path_to_raisim_activation_key = argv[1];
  idocp::QuadrupedSimulator<idocp::OCPSolver> simulator(path_to_raisim_activation_key, 
                                                        path_to_urdf_for_raisim, 
                                                        "../sim_result", "anymal");
  constexpr bool visualization = true;
  constexpr bool video_recording = false;
  simulator.run<MPCCallbackTrotting>(mpc, 1, 0.0025, 0, q, v, visualization, video_recording);
  // simulator.run<MPCCallbackStanding>(mpc, 0.75, 0.0025, 0, q, v, visualization, video_recording);
  mpc.getSolverHandle()->printSolution("end-effector", {14, 24, 34, 44});
  return 0;
}