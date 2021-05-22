#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_direction.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SplitDirectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ContactStatus& contact_status, 
                   const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, 
                           const ContactStatus& contact_status, 
                           const ImpulseStatus& impulse_status);

  double dt;
};


void SplitDirectionTest::test(const Robot& robot, const ContactStatus& contact_status,
                              const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  const int dimi = impulse_status.dimf();
  SplitDirection d(robot);
  EXPECT_EQ(d.dx.size(), dimx);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.du.size(), dimu);
  EXPECT_EQ(d.daf().size(), dimv);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dlmdgmm.size(), dimx);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.dnu_passive.size(), dim_passive);
  EXPECT_EQ(d.dxi().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimi(), 0);
  d.setContactStatus(contact_status);
  d.setImpulseStatus(impulse_status);
  EXPECT_EQ(d.dx.size(), dimx);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.du.size(), dimu);
  EXPECT_EQ(d.daf().size(), dimv+dimf);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.df().size(), dimf);
  EXPECT_EQ(d.dlmdgmm.size(), dimx);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dmu().size(), dimf);
  EXPECT_EQ(d.dnu_passive.size(), dim_passive);
  EXPECT_EQ(d.dxi().size(), dimi);
  EXPECT_EQ(d.dimf(), dimf);
  EXPECT_EQ(d.dimi(), dimi);
  d.setRandom();
  EXPECT_FALSE(d.dx.isZero());
  EXPECT_FALSE(d.daf().isZero());
  EXPECT_FALSE(d.dlmdgmm.isZero());
  EXPECT_FALSE(d.dbetamu().isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(d.dnu_passive.isZero());
  }
  if (dimi > 0) {
    EXPECT_FALSE(d.dxi().isZero());
  }
  EXPECT_TRUE(d.dx.head(dimv).isApprox(d.dq()));
  EXPECT_TRUE(d.dx.tail(dimv).isApprox(d.dv()));
  EXPECT_TRUE(d.daf().head(dimv).isApprox(d.da()));
  EXPECT_TRUE(d.daf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dlmdgmm.head(dimv).isApprox(d.dlmd()));
  EXPECT_TRUE(d.dlmdgmm.tail(dimv).isApprox(d.dgmm()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.dx.isZero());
  EXPECT_TRUE(d.daf().isZero());
  EXPECT_TRUE(d.dlmdgmm.isZero());
  EXPECT_TRUE(d.dbetamu().isZero());
  EXPECT_TRUE(d.dnu_passive.isZero());
  EXPECT_TRUE(d.dxi().isZero());
}


void SplitDirectionTest::testIsApprox(const Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  const int dimi = impulse_status.dimf();
  auto d = SplitDirection::Random(robot, contact_status, impulse_status);
  EXPECT_FALSE(d.dx.isZero());
  EXPECT_FALSE(d.daf().isZero());
  EXPECT_FALSE(d.dlmdgmm.isZero());
  EXPECT_FALSE(d.dbetamu().isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(d.dnu_passive.isZero());
  }
  if (dimi > 0) {
    EXPECT_FALSE(d.dxi().isZero());
  }
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
  d_ref.da().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.da() = d.da();
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
  if (contact_status.hasActiveContacts()) {
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
  if (impulse_status.hasActiveImpulse()) {
    d_ref.dxi().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.dxi() = d.dxi();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
  else {
    d_ref.dxi().setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
}


TEST_F(SplitDirectionTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.deactivateContact(0);
  impulse_status.activateImpulse(0);
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  impulse_status.activateImpulse(0);
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
}


TEST_F(SplitDirectionTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.deactivateContact(i);
  }
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}