#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/collision_checker.hpp"
#include "idocp/robot/robot.hpp"


namespace idocp {

class CollisionCheckerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void test(Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
};


void CollisionCheckerTest::test(Robot& robot) {
  CollisionChecker checker(robot);
  Eigen::VectorXd q(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  const auto is_collision = checker.check(robot, q);
  std::vector<Eigen::Vector3d> ee(robot.maxPointContacts(), Eigen::Vector3d::Zero());
  std::vector<bool> is_collision_ref(robot.maxPointContacts());
  robot.updateFrameKinematics(q);
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    ee[i] = robot.framePosition(robot.contactFramesIndices()[i]);
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
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(CollisionCheckerTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}