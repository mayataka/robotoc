#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseSplitDirectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);

};


void ImpulseSplitDirectionTest::test(const Robot& robot, const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimi = impulse_status.dimf();
  ImpulseSplitDirection d(robot);
  EXPECT_EQ(d.dx.size(), dimx);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.ddvf().size(), dimv);
  EXPECT_EQ(d.ddv().size(), dimv);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dlmdgmm.size(), dimx);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.dimi(), 0);
  d.setImpulseStatus(impulse_status);
  EXPECT_EQ(d.dx.size(), dimx);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.ddvf().size(), dimv+dimi);
  EXPECT_EQ(d.ddv().size(), dimv);
  EXPECT_EQ(d.df().size(), dimi);
  EXPECT_EQ(d.dlmdgmm.size(), dimx);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv+dimi);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dmu().size(), dimi);
  EXPECT_EQ(d.dimi(), dimi);
  d.setRandom();
  EXPECT_FALSE(d.dx.isZero());
  EXPECT_FALSE(d.ddvf().isZero());
  EXPECT_FALSE(d.dlmdgmm.isZero());
  EXPECT_FALSE(d.dbetamu().isZero());
  EXPECT_TRUE(d.dx.head(dimv).isApprox(d.dq()));
  EXPECT_TRUE(d.dx.tail(dimv).isApprox(d.dv()));
  EXPECT_TRUE(d.ddvf().head(dimv).isApprox(d.ddv()));
  EXPECT_TRUE(d.ddvf().tail(dimi).isApprox(d.df()));
  EXPECT_TRUE(d.dlmdgmm.head(dimv).isApprox(d.dlmd()));
  EXPECT_TRUE(d.dlmdgmm.tail(dimv).isApprox(d.dgmm()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimi).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.dx.isZero());
  EXPECT_TRUE(d.ddvf().isZero());
  EXPECT_TRUE(d.dlmdgmm.isZero());
  EXPECT_TRUE(d.dbetamu().isZero());
}


void ImpulseSplitDirectionTest::testIsApprox(const Robot& robot, 
                                             const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimi = impulse_status.dimf();
  auto d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dq().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dq() = d.dq();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dv().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dv() = d.dv();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.ddv().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.ddv() = d.ddv();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dlmd().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dlmd() = d.dlmd();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dgmm().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dgmm() = d.dgmm();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dbeta().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dbeta() = d.dbeta();
  EXPECT_TRUE(d.isApprox(d_ref));
  if (impulse_status.hasActiveImpulse()) {
    d_ref.df().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.df() = d.df();
    EXPECT_TRUE(d.isApprox(d_ref));
    d_ref.dmu().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.dmu() = d.dmu();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
  else {
    d_ref.df().setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
    d_ref.dmu().setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
}


TEST_F(ImpulseSplitDirectionTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(ImpulseSplitDirectionTest, floatingBase) {
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