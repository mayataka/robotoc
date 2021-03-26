#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/collision_checker.hpp"
#include "idocp/robot/robot.hpp"

#include "robot_factory.hpp"


namespace idocp {

class CollisionCheckerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(Robot& robot);
};


void CollisionCheckerTest::test(Robot& robot) {
  CollisionChecker checker(robot);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const auto is_collision = checker.check(robot, q);
  std::vector<Eigen::Vector3d> ee(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  std::vector<bool> is_collision_ref(robot.maxPointContacts());
  robot.updateFrameKinematics(q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    ee[i] = robot.framePosition(robot.contactFrames()[i]);
    if (ee[i].coeff(2) <= 0) is_collision_ref[i] = true;
    else is_collision_ref[i] = false;
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(is_collision[i], is_collision_ref[i]);
    EXPECT_TRUE(ee[i].isApprox(checker.contactFramePositions()[i]));
  }
  checker.printContactFramePositions();
}


TEST_F(CollisionCheckerTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(CollisionCheckerTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}