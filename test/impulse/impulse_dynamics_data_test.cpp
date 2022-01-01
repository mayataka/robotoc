#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_dynamics_data.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ImpulseDynamicsDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseDynamicsDataTest::test(const Robot& robot, 
                                   const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimi = impulse_status.dimi();
  ImpulseDynamicsData data(robot);
  data.setImpulseStatus(impulse_status);
  EXPECT_EQ(data.dImDddv.rows(), dimv);
  EXPECT_EQ(data.dImDddv.cols(), dimv);
  EXPECT_EQ(data.dImDCdqv().rows(), dimv+dimi);
  EXPECT_EQ(data.dImDCdqv().cols(), dimx);
  EXPECT_EQ(data.dImDCdq().rows(), dimv+dimi);
  EXPECT_EQ(data.dImDCdq().cols(), dimv);
  EXPECT_EQ(data.dCdq().rows(), dimi);
  EXPECT_EQ(data.dCdq().cols(), dimv);
  EXPECT_EQ(data.dCdv().rows(), dimi);
  EXPECT_EQ(data.dCdv().cols(), dimv);
  EXPECT_EQ(data.MJtJinv().rows(), dimv+dimi);
  EXPECT_EQ(data.MJtJinv().cols(), dimv+dimi);
  EXPECT_EQ(data.MJtJinv_dImDCdqv().rows(), dimv+dimi);
  EXPECT_EQ(data.MJtJinv_dImDCdqv().cols(), dimx);
  EXPECT_EQ(data.Qdvfqv().rows(), dimv+dimi);
  EXPECT_EQ(data.Qdvfqv().cols(), dimx);
  EXPECT_EQ(data.ImDC().size(), dimv+dimi);
  EXPECT_EQ(data.ImD().size(), dimv);
  EXPECT_EQ(data.C().size(), dimi);
  EXPECT_EQ(data.MJtJinv_ImDC().size(), dimv+dimi);
  EXPECT_EQ(data.ldvf().size(), dimv+dimi);
  EXPECT_EQ(data.ldv().size(), dimv);
  EXPECT_EQ(data.lf().size(), dimi);
  const Eigen::MatrixXd dImDddv_ref = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd dImDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimi, dimx);
  const Eigen::MatrixXd MJtJinv_ref = Eigen::MatrixXd::Random(dimv+dimi, dimv+dimi);
  const Eigen::MatrixXd MJtJinv_dImDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimi, dimx);
  const Eigen::MatrixXd Qdvfqv_ref = Eigen::MatrixXd::Random(dimv+dimi, dimx);
  const Eigen::VectorXd ImDC_ref = Eigen::VectorXd::Random(dimv+dimi);
  const Eigen::VectorXd MJtJinv_ImDC_ref = Eigen::VectorXd::Random(dimv+dimi);
  const Eigen::VectorXd ldvf_ref = Eigen::VectorXd::Random(dimv+dimi);
  data.dImDddv = dImDddv_ref;
  data.dImDCdqv() = dImDCdqv_ref;
  data.MJtJinv() = MJtJinv_ref;
  data.MJtJinv_dImDCdqv() = MJtJinv_dImDCdqv_ref;
  data.Qdvfqv() = Qdvfqv_ref;
  data.ImDC() = ImDC_ref;
  data.MJtJinv_ImDC() = MJtJinv_ImDC_ref;
  data.ldvf() = ldvf_ref;
  EXPECT_TRUE(data.dImDddv.isApprox(dImDddv_ref));
  EXPECT_TRUE(data.dImDCdqv().isApprox(dImDCdqv_ref));
  EXPECT_TRUE(data.dImDCdq().isApprox(dImDCdqv_ref.leftCols(dimv)));
  EXPECT_TRUE(data.dImDdq().isApprox(dImDCdqv_ref.block(0, 0, dimv, dimv)));
  EXPECT_TRUE(data.dCdq().isApprox(dImDCdqv_ref.block(dimv, 0, dimi, dimv)));
  EXPECT_TRUE(data.dCdv().isApprox(dImDCdqv_ref.block(dimv, dimv, dimi, dimv)));
  EXPECT_TRUE(data.MJtJinv().isApprox(MJtJinv_ref));
  EXPECT_TRUE(data.MJtJinv_dImDCdqv().isApprox(MJtJinv_dImDCdqv_ref));
  EXPECT_TRUE(data.Qdvfqv().isApprox(Qdvfqv_ref));
  EXPECT_TRUE(data.ImDC().isApprox(ImDC_ref));
  EXPECT_TRUE(data.ImD().isApprox(ImDC_ref.head(dimv)));
  EXPECT_TRUE(data.C().isApprox(ImDC_ref.tail(dimi)));
  EXPECT_TRUE(data.MJtJinv_ImDC().isApprox(MJtJinv_ImDC_ref));
  EXPECT_TRUE(data.ldvf().isApprox(ldvf_ref));
  EXPECT_TRUE(data.ldv().isApprox(ldvf_ref.head(dimv)));
  EXPECT_TRUE(data.lf().isApprox(ldvf_ref.tail(dimi)));
}


TEST_F(ImpulseDynamicsDataTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
}


TEST_F(ImpulseDynamicsDataTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}