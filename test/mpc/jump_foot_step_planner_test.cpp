#include <vector>

#include <gtest/gtest.h>

#include "robotoc/mpc/jump_foot_step_planner.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class JumpFootStepPlannerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(JumpFootStepPlannerTest, enableStancePhase) {
  const double dt = 0.05;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto planner = std::make_shared<JumpFootStepPlanner>(robot);
  const Eigen::Vector3d jump_length = Eigen::Vector3d::Random();
  const double jump_yaw = 0.1;
  planner->setJumpPattern(jump_length, jump_yaw);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto contact_status = robot.createContactStatus();
  contact_status.activateContact(0);
  contact_status.activateContact(1);
  contact_status.activateContact(2);
  contact_status.activateContact(3);
  const double t0 = 1.5;
  const int planning_steps = 10;
  planner->init(q);
  bool success = planner->plan(t0, q, v, contact_status, planning_steps);
  EXPECT_TRUE(success);
  const double t1 = t0 + dt;
  contact_status.deactivateContact(0);
  contact_status.deactivateContact(1);
  contact_status.deactivateContact(2);
  contact_status.deactivateContact(3);
  success = planner->plan(t1, q, v, contact_status, planning_steps);
  EXPECT_TRUE(success);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}