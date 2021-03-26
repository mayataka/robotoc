#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_backward_correction_data.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseSplitBackwardCorrectionDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((signed int) time(0));
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;
};


void ImpulseSplitBackwardCorrectionDataTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  int dimf = 3;
  if (robot.hasFloatingBase()) {
    std::random_device rnd;
    dimf += 3 * (rnd() % 4);
  }
  const int dimKKT = 4*robot.dimv() + 2*dimf;

  ImpulseSplitBackwardCorrectionData data(robot);
  data.setImpulseStatus(dimf);
  EXPECT_EQ(data.KKT_mat_inv().rows(), dimKKT);
  EXPECT_EQ(data.KKT_mat_inv().cols(), dimKKT);
  data.KKT_mat_inv().setRandom();
  EXPECT_TRUE(data.KKT_mat_inv().topLeftCorner(dimx, dimx).isApprox(data.auxMat()));
  EXPECT_EQ(data.splitDirection().size(), dimKKT);
  EXPECT_EQ(data.dlmd().size(), dimv);
  EXPECT_EQ(data.dgmm().size(), dimv);
  EXPECT_EQ(data.dmu().size(), dimf);
  EXPECT_EQ(data.df().size(), dimf);
  EXPECT_EQ(data.dq().size(), dimv);
  EXPECT_EQ(data.dv().size(), dimv);
  data.splitDirection().setRandom();
  EXPECT_TRUE(data.splitDirection().head(dimv).isApprox(data.dlmd()));
  EXPECT_TRUE(data.splitDirection().segment(dimv, dimv).isApprox(data.dgmm()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv, dimf).isApprox(data.dmu()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv, dimf).isApprox(data.dmu()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv+dimf, dimf).isApprox(data.df()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv+2*dimf, dimv).isApprox(data.dq()));
  EXPECT_TRUE(data.splitDirection().segment(3*dimv+2*dimf, dimv).isApprox(data.dv()));
}


TEST_F(ImpulseSplitBackwardCorrectionDataTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(ImpulseSplitBackwardCorrectionDataTest, floating_base) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}