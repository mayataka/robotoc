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

#include "idocp/utils/quadruped_simulator.hpp"


class MPCCallbackTrotting {
public:
  using FootStepTrottingCostPtr = std::shared_ptr<idocp::FootStepTrottingCost>;
  MPCCallbackTrotting(const idocp::Robot& robot, 
                      const FootStepTrottingCostPtr& LF_foot_cost, 
                      const FootStepTrottingCostPtr& LH_foot_cost, 
                      const FootStepTrottingCostPtr& RF_foot_cost, 
                      const FootStepTrottingCostPtr& RH_foot_cost)
    : robot_(robot),
      LF_foot_cost_(LF_foot_cost), 
      LH_foot_cost_(LH_foot_cost), 
      RF_foot_cost_(RF_foot_cost), 
      RH_foot_cost_(RH_foot_cost), 
      contact_points_even_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
      contact_points_odd_(robot.maxPointContacts(), Eigen::Vector3d::Zero()),
      contact_status_even_(robot.createContactStatus()),
      contact_status_odd_(robot.createContactStatus()),
      steps_(0),
      step_times_(),
      future_step_times_() {
    contact_status_even_.activateContacts({1, 2});
    contact_status_odd_.activateContacts({0, 3});
    for (int i=0; i<20; ++i) {
      future_step_times_.push_back(kStart+i*kPeriod);
    }
  } 

  template <typename OCPSolverType>
  void init(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            idocp::MPC<OCPSolverType>& mpc) {
    robot_.updateFrameKinematics(q);
    robot_.getContactPoints(contact_points_even_);
    mpc.getSolverHandle()->setContactPoints(0, contact_points_even_);
    contact_points_even_[0].coeffRef(2) = 0; // LF
    contact_points_even_[1].coeffRef(2) = 0; // LH
    contact_points_even_[2].coeffRef(2) = 0; // RF
    contact_points_even_[3].coeffRef(2) = 0; // RH
    LF_foot_cost_->set_q_3d_ref(contact_points_even_[0], kStepLength, kStepHeightFront);
    LH_foot_cost_->set_q_3d_ref(contact_points_even_[1], kStepLength, kStepHeightHip);
    RF_foot_cost_->set_q_3d_ref(contact_points_even_[2], kStepLength, kStepHeightFront);
    RH_foot_cost_->set_q_3d_ref(contact_points_even_[3], kStepLength, kStepHeightHip);
    LF_foot_cost_->set_period(kStart, kPeriod);
    LH_foot_cost_->set_period(kStart+kPeriod, kPeriod);
    RF_foot_cost_->set_period(kStart+kPeriod, kPeriod);
    RH_foot_cost_->set_period(kStart, kPeriod);
  }

  template <typename OCPSolverType>
  void callback(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v, idocp::MPC<OCPSolverType>& mpc) {
    if (!step_times_.empty() && step_times_.front() < t) {
      std::cout << "popFront()!!!" << std::endl;
      mpc.getSolverHandle()->popFrontDiscreteEvent();
      step_times_.pop_front();
      ++steps_;
    }
    if (steps_ % 2 == 0) {
      robot_.updateFrameKinematics(q);
      robot_.getContactPoints(contact_points_even_);
      mpc.getSolverHandle()->setContactPoints(0, contact_points_even_);
      if (step_times_.size() == 2) {
        std::cout << "update for odd" << std::endl;
        robot_.getContactPoints(contact_points_odd_);
        contact_points_odd_[1].coeffRef(0) += kStepLength;
        contact_points_odd_[2].coeffRef(0) += kStepLength;
        mpc.getSolverHandle()->setContactPoints(1, contact_points_odd_);
      }
      else if (future_step_times_.front() < t+kT) {
        std::cout << "pushBack for even!!!" << std::endl;
        robot_.getContactPoints(contact_points_even_);
        if (steps_ > 0) {
          contact_points_even_[1].coeffRef(0) += kStepLength;
          contact_points_even_[2].coeffRef(0) += kStepLength;
        }
        contact_status_even_.setContactPoints(contact_points_even_);
        mpc.getSolverHandle()->pushBackContactStatus(
            contact_status_even_, future_step_times_.front(), t);
        step_times_.push_back(future_step_times_.front());
        future_step_times_.pop_front();
      }
    }
    else {
      robot_.updateFrameKinematics(q);
      robot_.getContactPoints(contact_points_odd_);
      mpc.getSolverHandle()->setContactPoints(0, contact_points_odd_);
      if(step_times_.size() == 2) {
        std::cout << "update for even" << std::endl;
        robot_.getContactPoints(contact_points_even_);
        contact_points_odd_[0].coeffRef(0) += kStepLength;
        contact_points_odd_[3].coeffRef(0) += kStepLength;
        mpc.getSolverHandle()->setContactPoints(1, contact_points_even_);
      }
      else if (future_step_times_.front() < t+kT) {
        std::cout << "pushBack for odd!!!" << std::endl;
        robot_.getContactPoints(contact_points_odd_);
        if (steps_ > 0) {
          contact_points_odd_[0].coeffRef(0) += kStepLength;
          contact_points_odd_[3].coeffRef(0) += kStepLength;
        }
        contact_status_odd_.setContactPoints(contact_points_odd_);
        mpc.getSolverHandle()->pushBackContactStatus(
            contact_status_odd_, future_step_times_.front(), t);
        step_times_.push_back(future_step_times_.front());
        future_step_times_.pop_front();
      }
    }
    // update Cost
    robot_.getContactPoints(contact_points_even_);
    mpc.getSolverHandle()->setContactPoints(0, contact_points_even_);
    for (auto& e : contact_points_even_) {
      e.coeffRef(0) += kStepLength;
    }
    for (int i=0; i<100; ++i) {
      mpc.updateSolution(t, q, v);
    }
    mpc.computeKKTResidual(t, q, v);
    std::cout << "time t = " << t << ", KKTError = " << mpc.KKTError() << std::endl;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  idocp::Robot robot_;
  FootStepTrottingCostPtr LF_foot_cost_, LH_foot_cost_, RF_foot_cost_, RH_foot_cost_;
  std::vector<Eigen::Vector3d> contact_points_even_, contact_points_odd_;
  idocp::ContactStatus contact_status_even_, contact_status_odd_;
  int steps_;
  std::deque<double> step_times_, future_step_times_;
  static constexpr double kPeriod = 0.5;
  static constexpr double kStepLength = 0.15;
  static constexpr double kStepHeightFront = 0.1;
  static constexpr double kStepHeightHip = 0.1;
  static constexpr double kT = 0.5;
  static constexpr double kStart = 1.0;
  static constexpr double kSamplingPeriod = 0.025;
};



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
  robot.normalizeConfiguration(q_ref);
  Eigen::VectorXd v_ref(Eigen::VectorXd::Zero(robot.dimv()));
  v_ref <<  0.15, 0, 0, 0, 0, 0, 
              0, 0, 0,
              0, 0, 0,
              0, 0, 0,
              0, 0, 0;
  Eigen::VectorXd q_weight(Eigen::VectorXd::Zero(robot.dimv()));
  q_weight << 10, 10, 10, 10, 10, 10, 
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1,
               0.1, 0.1, 0.1;
  Eigen::VectorXd v_weight(Eigen::VectorXd::Zero(robot.dimv()));
  v_weight << 1, 1, 1, 1, 1, 1, 
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1,
              0.1, 0.1, 0.1;
  Eigen::VectorXd a_weight(Eigen::VectorXd::Zero(robot.dimv()));
  a_weight << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01,
              0.01, 0.01, 0.01;
  auto configuration_cost = std::make_shared<idocp::TimeVaryingConfigurationCost>(robot);
  const double t_ref_start = 1;
  Eigen::VectorXd q_test(Eigen::VectorXd::Zero(robot.dimq()));
  robot.integrateConfiguration(q_ref, v_ref, t_ref_start, q_test);
  configuration_cost->set_ref(t_ref_start, q_ref, v_ref);
  configuration_cost->set_q_weight(q_weight);
  configuration_cost->set_qf_weight(q_weight);
  configuration_cost->set_v_weight(v_weight);
  configuration_cost->set_vf_weight(v_weight);
  configuration_cost->set_a_weight(a_weight);
  cost->push_back(configuration_cost);

  auto LF_foot_cost = std::make_shared<idocp::FootStepTrottingCost>(robot, 14);
  auto LH_foot_cost = std::make_shared<idocp::FootStepTrottingCost>(robot, 24);
  auto RF_foot_cost = std::make_shared<idocp::FootStepTrottingCost>(robot, 34);
  auto RH_foot_cost = std::make_shared<idocp::FootStepTrottingCost>(robot, 44);
  const double weight = 10;
  LF_foot_cost->set_q_3d_weight(Eigen::Vector3d::Constant(weight));
  LF_foot_cost->set_qf_3d_weight(Eigen::Vector3d::Constant(weight));
  LH_foot_cost->set_q_3d_weight(Eigen::Vector3d::Constant(weight));
  LH_foot_cost->set_qf_3d_weight(Eigen::Vector3d::Constant(weight));
  RF_foot_cost->set_q_3d_weight(Eigen::Vector3d::Constant(weight));
  RF_foot_cost->set_qf_3d_weight(Eigen::Vector3d::Constant(weight));
  RH_foot_cost->set_q_3d_weight(Eigen::Vector3d::Constant(weight));
  RH_foot_cost->set_qf_3d_weight(Eigen::Vector3d::Constant(weight));
  cost->push_back(LF_foot_cost);
  cost->push_back(LH_foot_cost);
  cost->push_back(RF_foot_cost);
  cost->push_back(RH_foot_cost);

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
  auto impulse_configuration_cost = std::make_shared<idocp::ImpulseTimeVaryingConfigurationCost>(robot);
  impulse_configuration_cost->set_q_weight(q_weight);
  impulse_configuration_cost->set_ref(0, q_ref, v_ref);
  impulse_configuration_cost->set_v_weight(v_weight);
  impulse_configuration_cost->set_dv_weight(a_weight);
  cost->push_back(impulse_configuration_cost);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  impulse_force_cost->set_f_weight(f_weight);
  impulse_force_cost->set_f_ref(f_ref);
  cost->push_back(impulse_force_cost);

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
  const double T = 0.5;
  const int N = 20;
  const int max_num_impulse_phase = 2;
  const int num_proc = 4;
  const double t = 0;
  idocp::MPC<idocp::OCPSolver> mpc(robot, cost, constraints, T, N, max_num_impulse_phase, num_proc);

  Eigen::VectorXd q(Eigen::VectorXd::Zero(robot.dimq()));
  q << 0, 0, 0.4792, 0, 0, 0, 1, 
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

  mpc.initializeSolution(t, q, v, 100);
  mpc.computeKKTResidual(t, q, v);
  std::cout << mpc.KKTError() << std::endl;
  const std::string path_to_raisim_activation_key = argv[1];
  const std::string path_to_urdf_for_raisim = "../anymal/anymal_for_raisim.urdf";
  idocp::QuadrupedSimulator<idocp::OCPSolver> simulator(path_to_raisim_activation_key, 
                                                        path_to_urdf_for_raisim, 
                                                        "../sim_result", "trotting");
  constexpr bool visualization = true;
  constexpr bool video_recording = false;
  MPCCallbackTrotting mpc_callback(robot, LF_foot_cost, LH_foot_cost, RF_foot_cost, RH_foot_cost);
  simulator.run(mpc, mpc_callback, 1.5, 0.0025, 0, q, v, visualization, video_recording);
  return 0;
}