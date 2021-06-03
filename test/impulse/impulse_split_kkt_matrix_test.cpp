#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseSplitKKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseSplitKKTMatrixTest::test(const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTMatrix kkt_mat(robot);
  kkt_mat.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimi = impulse_status.dimf();
  const int dimx = 2 * robot.dimv();
  EXPECT_EQ(kkt_mat.dimi(), dimi);
  EXPECT_EQ(kkt_mat.Fxx.rows(), 2*dimv);
  EXPECT_EQ(kkt_mat.Fxx.cols(), 2*dimv);
  EXPECT_EQ(kkt_mat.Fqq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fqq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fqv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fqv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fvq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fvq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Fvv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Fvv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qxx.rows(), 2*dimv);
  EXPECT_EQ(kkt_mat.Qxx.cols(), 2*dimv);
  EXPECT_EQ(kkt_mat.Qqq().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qqv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qvq().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qvv().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qvv().cols(), dimv);
  EXPECT_EQ(kkt_mat.Qdvdv.rows(), dimv);
  EXPECT_EQ(kkt_mat.Qdvdv.cols(), dimv);
  EXPECT_EQ(kkt_mat.Qff().rows(), dimi);
  EXPECT_EQ(kkt_mat.Qff().cols(), dimi);
  EXPECT_EQ(kkt_mat.Qqf().rows(), dimv);
  EXPECT_EQ(kkt_mat.Qqf().cols(), dimi);
  if (robot.hasFloatingBase()) {
    EXPECT_EQ(kkt_mat.Fqq_prev.rows(), dimv);
    EXPECT_EQ(kkt_mat.Fqq_prev.cols(), dimv);
  }
  else {
    EXPECT_EQ(kkt_mat.Fqq_prev.rows(), 0);
    EXPECT_EQ(kkt_mat.Fqq_prev.cols(), 0);
  }
  EXPECT_TRUE(kkt_mat.Fxx.isZero());
  EXPECT_TRUE(kkt_mat.Qxx.isZero());
  EXPECT_TRUE(kkt_mat.Qdvdv.isZero());
  EXPECT_TRUE(kkt_mat.Qff().isZero());
  EXPECT_TRUE(kkt_mat.Qqf().isZero());
  EXPECT_TRUE(kkt_mat.Fqq_prev.isZero());

  kkt_mat.Fxx.setRandom();
  kkt_mat.Qxx.setRandom();
  EXPECT_TRUE(kkt_mat.Fxx.topLeftCorner(dimv, dimv).isApprox(kkt_mat.Fqq()));
  EXPECT_TRUE(kkt_mat.Fxx.topRightCorner(dimv, dimv).isApprox(kkt_mat.Fqv()));
  EXPECT_TRUE(kkt_mat.Fxx.bottomLeftCorner(dimv, dimv).isApprox(kkt_mat.Fvq()));
  EXPECT_TRUE(kkt_mat.Fxx.bottomRightCorner(dimv, dimv).isApprox(kkt_mat.Fvv()));
  EXPECT_TRUE(kkt_mat.Qxx.topLeftCorner(dimv, dimv).isApprox(kkt_mat.Qqq()));
  EXPECT_TRUE(kkt_mat.Qxx.topRightCorner(dimv, dimv).isApprox(kkt_mat.Qqv()));
  EXPECT_TRUE(kkt_mat.Qxx.bottomLeftCorner(dimv, dimv).isApprox(kkt_mat.Qvq()));
  EXPECT_TRUE(kkt_mat.Qxx.bottomRightCorner(dimv, dimv).isApprox(kkt_mat.Qvv()));
}


void ImpulseSplitKKTMatrixTest::testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTMatrix kkt_mat(robot);
  kkt_mat.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimi = impulse_status.dimf();
  const int dimx = 2 * robot.dimv();
  kkt_mat.Fxx.setRandom();
  kkt_mat.Qxx.setRandom();
  kkt_mat.Qdvdv.setRandom();
  kkt_mat.Qff().setRandom();
  kkt_mat.Qqf().setRandom();
  auto kkt_mat_ref = kkt_mat;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fxx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Fxx = kkt_mat.Fxx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxx.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qxx = kkt_mat.Qxx;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qdvdv.setRandom();
  EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
  kkt_mat_ref.Qdvdv = kkt_mat.Qdvdv;
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  if (impulse_status.hasActiveImpulse()) {
    kkt_mat_ref.Qff().setRandom();
    EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qff() = kkt_mat.Qff();
    EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qqf().setRandom();
    EXPECT_FALSE(kkt_mat.isApprox(kkt_mat_ref));
    kkt_mat_ref.Qqf() = kkt_mat.Qqf();
    EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  }
  else {
    kkt_mat_ref.Qff().setRandom();
    kkt_mat_ref.Qqf().setRandom();
    EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
  }
}


TEST_F(ImpulseSplitKKTMatrixTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(ImpulseSplitKKTMatrixTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}