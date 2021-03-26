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

  static void testSize(const Robot& robot, 
                       const ContactStatus& contact_status, 
                       const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, 
                           const ContactStatus& contact_status, 
                           const ImpulseStatus& impulse_status);

  double dt;
};


void SplitDirectionTest::testSize(const Robot& robot, 
                                  const ContactStatus& contact_status,
                                  const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dimi = impulse_status.dimf();
  SplitDirection d(robot);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dxi().size(), 0);
  EXPECT_EQ(d.du().size(), dimu);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dnu_passive.size(), 6);
  EXPECT_EQ(d.daf().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  d.setContactStatus(contact_status);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.du().size(), dimu);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.dxi().size(), 0);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.dnu_passive.size(), 6);
  EXPECT_EQ(d.daf().size(), dimv+dimf);
  EXPECT_EQ(d.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d.df().size(), dimf);
  EXPECT_EQ(d.dmu().size(), dimf);
  EXPECT_EQ(d.dimf(), dimf);
  d.setImpulseStatus(impulse_status);
  EXPECT_EQ(d.dxi().size(), dimi);
  EXPECT_EQ(d.dimi(), dimi);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.split_direction.size());
  d.split_direction = split_direction;
  const Eigen::VectorXd dlmd = split_direction.segment(0,  dimv);
  const Eigen::VectorXd dgmm = split_direction.segment(dimv,  dimv);
  const Eigen::VectorXd du   = split_direction.segment(2*dimv,  dimu);
  const Eigen::VectorXd dq   = split_direction.segment(2*dimv+dimu,  dimv);
  const Eigen::VectorXd dv   = split_direction.segment(2*dimv+dimu+dimv,  dimv);
  const Eigen::VectorXd dx   = split_direction.segment(2*dimv+dimu, 2*dimv);
  const Eigen::VectorXd dxi  = Eigen::VectorXd::Random(dimi);
  d.dxi() = dxi;
  const Eigen::VectorXd da   = Eigen::VectorXd::Random(dimv);
  d.da() = da;
  const Eigen::VectorXd df = Eigen::VectorXd::Random(dimf);
  d.df() = df;
  const Eigen::VectorXd dmu = Eigen::VectorXd::Random(dimf);
  d.dmu() = dmu;
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(dimv);
  d.dbeta() = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(du.isApprox(d.du()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(dxi.isApprox(d.dxi()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(d.daf().head(dimv).isApprox(d.da()));
  EXPECT_TRUE(d.daf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.split_direction.isZero());
  EXPECT_TRUE(d.dlmd().isZero());
  EXPECT_TRUE(d.dgmm().isZero());
  EXPECT_TRUE(d.du().isZero());
  EXPECT_TRUE(d.dq().isZero());
  EXPECT_TRUE(d.dv().isZero());
  EXPECT_TRUE(d.dx().isZero());
  EXPECT_TRUE(d.dxi().isZero());
  EXPECT_TRUE(d.da().isZero());
  EXPECT_TRUE(d.df().isZero());
  EXPECT_TRUE(d.daf().isZero());
  EXPECT_TRUE(d.dbeta().isZero());
  EXPECT_TRUE(d.dmu().isZero());
  EXPECT_TRUE(d.dbetamu().isZero());
  EXPECT_TRUE(d.dnu_passive.isZero());
  d.setRandom();
  EXPECT_FALSE(d.split_direction.isZero());
  if (dimi > 0) {
    EXPECT_FALSE(d.dxi().isZero());
  }
  EXPECT_FALSE(d.dlmd().isZero());
  EXPECT_FALSE(d.dgmm().isZero());
  EXPECT_FALSE(d.du().isZero());
  EXPECT_FALSE(d.dq().isZero());
  EXPECT_FALSE(d.dv().isZero());
  EXPECT_FALSE(d.dx().isZero());
  EXPECT_FALSE(d.da().isZero());
  if (dimf > 0) {
    EXPECT_FALSE(d.df().isZero());
  }
  EXPECT_FALSE(d.daf().isZero());
  EXPECT_FALSE(d.dbeta().isZero());
  if (dimf > 0) {
    EXPECT_FALSE(d.dmu().isZero());
  }
  EXPECT_FALSE(d.dbetamu().isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(d.dnu_passive.isZero());
  }
  EXPECT_TRUE(d.daf().head(dimv).isApprox(d.da()));
  EXPECT_TRUE(d.daf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  const SplitDirection d_random = SplitDirection::Random(robot, contact_status, 
                                                         impulse_status);
  EXPECT_EQ(d_random.dlmd().size(), dimv);
  EXPECT_EQ(d_random.dgmm().size(), dimv);
  EXPECT_EQ(d_random.du().size(), dimu);
  EXPECT_EQ(d_random.dq().size(), dimv);
  EXPECT_EQ(d_random.dv().size(), dimv);
  EXPECT_EQ(d_random.dx().size(), dimx);
  EXPECT_EQ(d_random.dxi().size(), dimi);
  EXPECT_EQ(d_random.da().size(), dimv);
  EXPECT_EQ(d_random.dbeta().size(), dimv);
  EXPECT_EQ(d_random.dnu_passive.size(), 6);
  EXPECT_EQ(d_random.daf().size(), dimv+dimf);
  EXPECT_EQ(d_random.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d_random.df().size(), dimf);
  EXPECT_EQ(d_random.dmu().size(), dimf);
  EXPECT_EQ(d_random.dimf(), dimf);
  EXPECT_EQ(d_random.dimi(), dimi);
  EXPECT_FALSE(d_random.split_direction.isZero());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.du().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.dx().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(d_random.dnu_passive.isZero());
  }
  if (dimf > 0) {
    EXPECT_FALSE(d_random.df().isZero());
  }
  EXPECT_FALSE(d_random.daf().isZero());
  EXPECT_FALSE(d_random.dbeta().isZero());
  if (dimf > 0) {
    EXPECT_FALSE(d_random.dmu().isZero());
  }
  if (dimi > 0) {
    EXPECT_FALSE(d_random.dxi().isZero());
  }
  EXPECT_FALSE(d_random.dbetamu().isZero());
  EXPECT_TRUE(d_random.daf().head(dimv).isApprox(d_random.da()));
  EXPECT_TRUE(d_random.daf().tail(dimf).isApprox(d_random.df()));
  EXPECT_TRUE(d_random.dbetamu().head(dimv).isApprox(d_random.dbeta()));
  EXPECT_TRUE(d_random.dbetamu().tail(dimf).isApprox(d_random.dmu()));
}


void SplitDirectionTest::testIsApprox(const Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dimi = impulse_status.dimf();
  SplitDirection d(robot);
  d.setRandom(contact_status);
  d.setImpulseStatus(impulse_status);
  d.dxi().setRandom();
  EXPECT_FALSE(d.dlmd().isZero());
  EXPECT_FALSE(d.dgmm().isZero());
  EXPECT_FALSE(d.du().isZero());
  EXPECT_FALSE(d.dq().isZero());
  EXPECT_FALSE(d.dv().isZero());
  EXPECT_FALSE(d.dx().isZero());
  EXPECT_FALSE(d.da().isZero());
  EXPECT_FALSE(d.dbeta().isZero());
  if (dimi > 0) {
    EXPECT_FALSE(d.dxi().isZero());
  }
  SplitDirection d_ref = d;
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dlmd().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dlmd() = d.dlmd();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dgmm().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dgmm() = d.dgmm();
  EXPECT_TRUE(d.isApprox(d_ref));
  if (dimi > 0) {
    d_ref.dxi().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.dxi() = d.dxi();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
  d_ref.du().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.du() = d.du();
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
  if (robot.hasFloatingBase()) {
    d_ref.dnu_passive.setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.dnu_passive = d.dnu_passive;
    EXPECT_TRUE(d.isApprox(d_ref));
  }
  else {
    d_ref.dnu_passive.setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
}


TEST_F(SplitDirectionTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.deactivateContact(0);
  impulse_status.activateImpulse(0);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  impulse_status.activateImpulse(0);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
}


TEST_F(SplitDirectionTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.deactivateContact(i);
  }
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}