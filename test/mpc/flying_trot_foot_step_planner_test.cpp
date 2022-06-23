#include <vector>

#include <gtest/gtest.h>

#include "robotoc/mpc/flying_trot_foot_step_planner.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class FlyingTrotFootStepPlannerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(FlyingTrotFootStepPlannerTest, test) {
  const double dt = 0.05;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto planner = std::make_shared<FlyingTrotFootStepPlanner>(robot);
  const Eigen::Vector3d vcom_cmd = Eigen::Vector3d::Random();
  const double yaw_rate_cmd = 0.01;
  const double swing_time = 0.5;
  const double stance_time = 0.1;
  const double gain = 0.7;
  planner->setRaibertGaitPattern(vcom_cmd, yaw_rate_cmd, swing_time, stance_time, gain);
  const Eigen::Vector3d step_length = 2 * (swing_time+stance_time) * vcom_cmd;
  const double step_yaw = swing_time * yaw_rate_cmd;
  planner->setGaitPattern(step_length, step_yaw);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto contact_status = robot.createContactStatus();
  contact_status.activateContact(0);
  contact_status.deactivateContact(1);
  contact_status.deactivateContact(2);
  contact_status.activateContact(3);
  const double t0 = 1.5;
  const int planning_steps = 10;
  planner->init(q);
  bool success = planner->plan(t0, q, v, contact_status, planning_steps);
  EXPECT_TRUE(success);
  const double t1 = t0 + dt;
  success = planner->plan(t1, q, v, contact_status, planning_steps);
  EXPECT_TRUE(success);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}