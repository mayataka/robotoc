#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_backward_correction_data.hpp"

#include "test_helper.hpp"

namespace idocp {

class SplitBackwardCorrectionDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((signed int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


void SplitBackwardCorrectionDataTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimKKT = 4*robot.dimv() + robot.dimu();

  SplitBackwardCorrectionData data(robot);
  EXPECT_EQ(data.KKT_mat_inv().rows(), dimKKT);
  EXPECT_EQ(data.KKT_mat_inv().cols(), dimKKT);
  data.KKT_mat_inv().setRandom();
  EXPECT_TRUE(data.KKT_mat_inv().topLeftCorner(dimx, dimx).isApprox(data.auxMat()));
  EXPECT_EQ(data.splitDirection().size(), dimKKT);
  EXPECT_EQ(data.dlmd().size(), dimv);
  EXPECT_EQ(data.dgmm().size(), dimv);
  EXPECT_EQ(data.dxi().size(), 0);
  EXPECT_EQ(data.du().size(), dimu);
  EXPECT_EQ(data.dq().size(), dimv);
  EXPECT_EQ(data.dv().size(), dimv);
  data.splitDirection().setRandom();
  EXPECT_TRUE(data.splitDirection().head(dimv).isApprox(data.dlmd()));
  EXPECT_TRUE(data.splitDirection().segment(dimv, dimv).isApprox(data.dgmm()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv, dimu).isApprox(data.du()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv+dimu, dimv).isApprox(data.dq()));
  EXPECT_TRUE(data.splitDirection().segment(3*dimv+dimu, dimv).isApprox(data.dv()));
  int dimi = 3;
  if (robot.hasFloatingBase()) {
    std::random_device rnd;
    dimi += 3 * (rnd() % 4);
  }
  data.setImpulseStatus(dimi);
  EXPECT_EQ(data.KKT_mat_inv().rows(), dimKKT+dimi);
  EXPECT_EQ(data.KKT_mat_inv().cols(), dimKKT+dimi);
  data.KKT_mat_inv().setRandom();
  EXPECT_TRUE(data.KKT_mat_inv().topLeftCorner(dimx, dimx).isApprox(data.auxMat()));
  EXPECT_EQ(data.splitDirection().size(), dimKKT+dimi);
  EXPECT_EQ(data.dlmd().size(), dimv);
  EXPECT_EQ(data.dgmm().size(), dimv);
  EXPECT_EQ(data.dxi().size(), dimi);
  EXPECT_EQ(data.du().size(), dimu);
  EXPECT_EQ(data.dq().size(), dimv);
  EXPECT_EQ(data.dv().size(), dimv);
  data.splitDirection().setRandom();
  EXPECT_TRUE(data.splitDirection().head(dimv).isApprox(data.dlmd()));
  EXPECT_TRUE(data.splitDirection().segment(dimv, dimv).isApprox(data.dgmm()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv, dimi).isApprox(data.dxi()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv+dimi, dimu).isApprox(data.du()));
  EXPECT_TRUE(data.splitDirection().segment(2*dimv+dimi+dimu, dimv).isApprox(data.dq()));
  EXPECT_TRUE(data.splitDirection().segment(3*dimv+dimi+dimu, dimv).isApprox(data.dv()));
}


TEST_F(SplitBackwardCorrectionDataTest, fixedBase) {
  test(fixed_base_robot);
}


TEST_F(SplitBackwardCorrectionDataTest, floating_base) {
  test(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}