#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"

#include "robot_factory.hpp"

namespace idocp {

class ImpulseSplitKKTResidualTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }
  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void test_isApprox(const Robot& robot, 
                           const ImpulseStatus& impulse_status);

  virtual void TearDown() {
  }
};


void ImpulseSplitKKTResidualTest::test(const Robot& robot, const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2 * robot.dimv();
  const int dimi = impulse_status.dimf();
  EXPECT_EQ(kkt_res.dimi(), dimi);
  EXPECT_EQ(kkt_res.Fx.size(), dimx);
  EXPECT_EQ(kkt_res.Fq().size(), dimv);
  EXPECT_EQ(kkt_res.Fv().size(), dimv);
  EXPECT_EQ(kkt_res.lx.size(), dimx);
  EXPECT_EQ(kkt_res.lq().size(), dimv);
  EXPECT_EQ(kkt_res.lv().size(), dimv);
  EXPECT_EQ(kkt_res.ldv.size(), dimv);
  EXPECT_EQ(kkt_res.lf().size(), dimi);
  EXPECT_EQ(kkt_res.dimi(), dimi);
  EXPECT_TRUE(kkt_res.Fx.isZero());
  EXPECT_TRUE(kkt_res.lx.isZero());
  EXPECT_TRUE(kkt_res.ldv.isZero());
  EXPECT_TRUE(kkt_res.lf().isZero());
  kkt_res.Fx.setRandom();
  kkt_res.lx.setRandom();

  EXPECT_TRUE(kkt_res.Fx.head(dimv).isApprox(kkt_res.Fq()));
  EXPECT_TRUE(kkt_res.Fx.tail(dimv).isApprox(kkt_res.Fv()));
  EXPECT_TRUE(kkt_res.lx.head(dimv).isApprox(kkt_res.lq()));
  EXPECT_TRUE(kkt_res.lx.tail(dimv).isApprox(kkt_res.lv()));

  kkt_res.ldv.setRandom();
  kkt_res.lf().setRandom();
  const double err = kkt_res.KKTError();
  const double err_ref = kkt_res.Fx.squaredNorm() 
                          + kkt_res.lx.squaredNorm()
                          + kkt_res.ldv.squaredNorm()
                          + kkt_res.lf().squaredNorm();
  EXPECT_DOUBLE_EQ(err, err_ref);

  const double vio = kkt_res.constraintViolation();
  const double vio_ref = kkt_res.Fx.template lpNorm<1>();
  EXPECT_DOUBLE_EQ(vio, vio_ref);
}


void ImpulseSplitKKTResidualTest::test_isApprox(const Robot& robot, 
                                               const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_res.Fx.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.ldv.setRandom();
  kkt_res.lf().setRandom();
  auto kkt_res_ref = kkt_res;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.Fx.setRandom();
  EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.Fx = kkt_res_ref.Fx;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.lx.setRandom();
  EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.lx = kkt_res_ref.lx;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.ldv.setRandom();
  EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.ldv = kkt_res_ref.ldv;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  if (impulse_status.hasActiveImpulse()) {
    kkt_res_ref.lf().setRandom();
    EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
    kkt_res.lf() = kkt_res_ref.lf();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  else {
    kkt_res_ref.lf().setRandom();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
}


TEST_F(ImpulseSplitKKTResidualTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  test_isApprox(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
  test_isApprox(robot, impulse_status);
}


TEST_F(ImpulseSplitKKTResidualTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  test_isApprox(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
  test_isApprox(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}