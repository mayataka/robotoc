#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class RiccatiFactorizationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot);
};


void RiccatiFactorizationTest::test(const Robot& robot) {
  ASSERT_TRUE(robot.maxPointContacts() == 0);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  SplitRiccatiFactorization riccati(robot);
  EXPECT_EQ(riccati.P.rows(), dimx);
  EXPECT_EQ(riccati.P.cols(), dimx);
  EXPECT_EQ(riccati.s.size(), dimx);
  EXPECT_EQ(riccati.Pqq().rows(), dimv);
  EXPECT_EQ(riccati.Pqq().cols(), dimv);
  EXPECT_EQ(riccati.Pqv().rows(), dimv);
  EXPECT_EQ(riccati.Pqv().cols(), dimv);
  EXPECT_EQ(riccati.Pvq().rows(), dimv);
  EXPECT_EQ(riccati.Pvq().cols(), dimv);
  EXPECT_EQ(riccati.Pvv().rows(), dimv);
  EXPECT_EQ(riccati.Pvv().cols(), dimv);
  EXPECT_EQ(riccati.sq().size(), dimv);
  EXPECT_EQ(riccati.sv().size(), dimv);
  EXPECT_EQ(riccati.Gmm.size(), dimx);
}


TEST_F(RiccatiFactorizationTest, fixed_base) {
  auto robot = testhelper::CreateFixedBaseRobot();
  test(robot);
}


TEST_F(RiccatiFactorizationTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot();
  test(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}