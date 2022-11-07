#include <vector>

#include <gtest/gtest.h>

#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/local_contact_force_cost.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/solver/solver_options.hpp"
#include "robotoc/mpc/control_policy.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ControlPolicyTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(ControlPolicyTest, test) {
  const double baumgarte_time_step = 0.5 / 20;
  auto robot = testhelper::CreateQuadrupedalRobot(baumgarte_time_step);
  const int LF_foot_id = 12;
  const int LH_foot_id = 22;
  const int RF_foot_id = 32;
  const int RH_foot_id = 42;
  const std::vector<int> contact_frames = {LF_foot_id, LH_foot_id, RF_foot_id, RH_foot_id}; 

  // Create a cost function.
  auto cost = std::make_shared<robotoc::CostFunction>();
  Eigen::VectorXd q_standing(robot.dimq());
  q_standing << 0, 0, 0.4792, 0, 0, 0, 1, 
                -0.1,  0.7, -1.0, 
                -0.1, -0.7,  1.0, 
                 0.1,  0.7, -1.0, 
                 0.1, -0.7,  1.0;
  Eigen::VectorXd v_ref(robot.dimv());
  v_ref << 0, 0, 0, 0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0, 
           0, 0, 0;
  auto config_cost = std::make_shared<robotoc::ConfigurationSpaceCost>(robot);
  config_cost->set_q_weight(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_q_ref(q_standing);
  config_cost->set_q_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 10));
  config_cost->set_v_weight(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_v_weight_terminal(Eigen::VectorXd::Constant(robot.dimv(), 1));
  config_cost->set_a_weight(Eigen::VectorXd::Constant(robot.dimv(), 0.01));
  cost->add("config_cost", config_cost);
  auto local_contact_force_cost = std::make_shared<robotoc::LocalContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight, f_ref;
  for (int i=0; i<contact_frames.size(); ++i) {
    Eigen::Vector3d fw; 
    fw << 0.001, 0.001, 0.001;
    f_weight.push_back(fw);
    Eigen::Vector3d fr; 
    fr << 0, 0, 70;
    f_ref.push_back(fr);
  }
  local_contact_force_cost->set_f_weight(f_weight);
  local_contact_force_cost->set_f_ref(f_ref);
  cost->add("local_contact_force_cost", local_contact_force_cost);

  // Create inequality constraints.
  auto constraints = std::make_shared<robotoc::Constraints>();
  auto joint_position_lower = std::make_shared<robotoc::JointPositionLowerLimit>(robot);
  auto joint_position_upper = std::make_shared<robotoc::JointPositionUpperLimit>(robot);
  auto joint_velocity_lower = std::make_shared<robotoc::JointVelocityLowerLimit>(robot);
  auto joint_velocity_upper = std::make_shared<robotoc::JointVelocityUpperLimit>(robot);
  auto joint_torques_lower  = std::make_shared<robotoc::JointTorquesLowerLimit>(robot);
  auto joint_torques_upper  = std::make_shared<robotoc::JointTorquesUpperLimit>(robot);
  auto friction_cone        = std::make_shared<robotoc::FrictionCone>(robot);
  constraints->add("joint_position_lower", joint_position_lower);
  constraints->add("joint_position_upper", joint_position_upper);
  constraints->add("joint_velocity_lower", joint_velocity_lower);
  constraints->add("joint_velocity_upper", joint_velocity_upper);
  constraints->add("joint_torques_lower", joint_torques_lower);
  constraints->add("joint_torques_upper", joint_torques_upper);
  constraints->add("friction_cone", friction_cone);

  // Create the contact sequence
  auto contact_sequence = std::make_shared<robotoc::ContactSequence>(robot);

  auto contact_status_standing = robot.createContactStatus();
  contact_status_standing.activateContacts({0, 1, 2, 3});
  robot.updateFrameKinematics(q_standing);
  const std::vector<Eigen::Vector3d> contact_positions = {robot.framePosition(LF_foot_id), 
                                                          robot.framePosition(LH_foot_id),
                                                          robot.framePosition(RF_foot_id),
                                                          robot.framePosition(RH_foot_id)};
  contact_status_standing.setContactPlacements(contact_positions);
  contact_sequence->init(contact_status_standing);

  // Create OCPSolver
  const double T = 0.5;
  const int N = 20;
  robotoc::OCP ocp(robot, cost, constraints, contact_sequence, T, N);
  auto solver_options = robotoc::SolverOptions();
  solver_options.nthreads = 4;
  robotoc::OCPSolver ocp_solver(ocp, solver_options);

  // Initial time and initial state
  const double t = 0;
  const Eigen::VectorXd q = q_standing;
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());

  ocp_solver.discretize(t);
  ocp_solver.setSolution("q", q);
  ocp_solver.setSolution("v", v);
  Eigen::Vector3d f_init;
  f_init << 0, 0, 0.25*robot.totalWeight();
  ocp_solver.setSolution("f", f_init);

  auto contact_status_flying = robot.createContactStatus();
  contact_sequence->push_back(contact_status_flying, 0.2);

  contact_sequence->push_back(contact_status_standing, 0.4);

  auto options = SolverOptions();
  options.max_iter = 10;
  ocp_solver.setSolverOptions(options);
  ocp_solver.solve(t, q, v);

  const double t0 = t - T / N * 5;
  const int nJ = robot.dimu();
  for (int i=0; i<100; ++i) {
    std::cout << i << std::endl;
    const auto policy = ControlPolicy(ocp_solver, t0+0.01);
    EXPECT_FALSE(policy.tauJ.hasNaN());
    EXPECT_FALSE(policy.qJ.hasNaN());
    EXPECT_FALSE(policy.dqJ.hasNaN());
    EXPECT_FALSE(policy.Kp.hasNaN());
    EXPECT_FALSE(policy.Kd.hasNaN());
    EXPECT_FALSE(policy.tauJ.isZero());
    EXPECT_FALSE(policy.qJ.isZero());
    EXPECT_FALSE(policy.dqJ.isZero());
    EXPECT_FALSE(policy.Kp.isZero());
    EXPECT_FALSE(policy.Kd.isZero());
    EXPECT_EQ(policy.tauJ.size(), nJ);
    EXPECT_EQ(policy.qJ.size(), nJ);
    EXPECT_EQ(policy.dqJ.size(), nJ);
    EXPECT_EQ(policy.Kp.rows(), nJ);
    EXPECT_EQ(policy.Kp.cols(), nJ);
    EXPECT_EQ(policy.Kd.rows(), nJ);
    EXPECT_EQ(policy.Kd.cols(), nJ);
  }

}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}