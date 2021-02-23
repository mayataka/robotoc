#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"


namespace idocp {

class RiccatiFactorizationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
    fixed_base_robot_contact = Robot(fixed_base_urdf, {18});
    floating_base_robot_contact = Robot(floating_base_urdf, {14, 24, 34, 44});
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
  Robot fixed_base_robot_contact, floating_base_robot_contact;
};


void RiccatiFactorizationTest::test(const Robot& robot) {
  ASSERT_TRUE(robot.maxPointContacts() == 0);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  SplitRiccatiFactorization riccati(robot);
  EXPECT_EQ(riccati.Pqq.rows(), dimv);
  EXPECT_EQ(riccati.Pqq.cols(), dimv);
  EXPECT_EQ(riccati.Pqv.rows(), dimv);
  EXPECT_EQ(riccati.Pqv.cols(), dimv);
  EXPECT_EQ(riccati.Pvq.rows(), dimv);
  EXPECT_EQ(riccati.Pvq.cols(), dimv);
  EXPECT_EQ(riccati.Pvv.rows(), dimv);
  EXPECT_EQ(riccati.Pvv.cols(), dimv);
  EXPECT_EQ(riccati.sq.size(), dimv);
  EXPECT_EQ(riccati.sv.size(), dimv);
}


TEST_F(RiccatiFactorizationTest, fixed_base) {
  test(fixed_base_robot);
}


TEST_F(RiccatiFactorizationTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}