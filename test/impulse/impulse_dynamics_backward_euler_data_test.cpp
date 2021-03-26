#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler_data.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseDynamicsForwardEulerDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseDynamicsForwardEulerDataTest::test(const Robot& robot, 
                                               const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimf();
  ImpulseDynamicsBackwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  EXPECT_EQ(data.dImDdq.rows(), dimv);
  EXPECT_EQ(data.dImDdq.cols(), dimv);
  EXPECT_EQ(data.dImDddv.rows(), dimv);
  EXPECT_EQ(data.dImDddv.cols(), dimv);
  EXPECT_EQ(data.Minv.rows(), dimv);
  EXPECT_EQ(data.Minv.cols(), dimv);
  EXPECT_EQ(data.Qdvq.rows(), dimv);
  EXPECT_EQ(data.Qdvq.cols(), dimv);
  EXPECT_EQ(data.ImD.size(), dimv);
  EXPECT_EQ(data.Minv_ImD.size(), dimv);
  EXPECT_EQ(data.ldv.size(), dimv);
  EXPECT_EQ(data.Qdvf().rows(), dimv);
  EXPECT_EQ(data.Qdvf().cols(), dimf);

  const Eigen::MatrixXd Qdvf_ref = Eigen::MatrixXd::Random(dimv, dimf);
  data.Qdvf() = Qdvf_ref;
  EXPECT_TRUE(data.Qdvf().isApprox(Qdvf_ref));
}


TEST_F(ImpulseDynamicsForwardEulerDataTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
}


TEST_F(ImpulseDynamicsForwardEulerDataTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}