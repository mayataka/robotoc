#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SplitKKTResidualTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  static void test(const Robot& robot, 
                       const ContactStatus& contact_status, 
                       const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, 
                           const ContactStatus& contact_status, 
                           const ImpulseStatus& impulse_status);

  virtual void TearDown() {
  }

  double dt;
};


void SplitKKTResidualTest::test(const Robot& robot, 
                                const ContactStatus& contact_status,
                                const ImpulseStatus& impulse_status) {
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dimi = impulse_status.dimf();
  const int dim_passive = robot.dim_passive();
  EXPECT_EQ(kkt_res.dimf(), dimf);
  EXPECT_EQ(kkt_res.dimi(), dimi);
  EXPECT_EQ(kkt_res.Fx.size(), dimx);
  EXPECT_EQ(kkt_res.Fq().size(), dimv);
  EXPECT_EQ(kkt_res.Fv().size(), dimv);
  EXPECT_EQ(kkt_res.P().size(), dimi);
  EXPECT_EQ(kkt_res.lx.size(), dimx);
  EXPECT_EQ(kkt_res.lq().size(), dimv);
  EXPECT_EQ(kkt_res.lv().size(), dimv);
  EXPECT_EQ(kkt_res.la.size(), dimv);
  EXPECT_EQ(kkt_res.lu.size(), dimu);
  EXPECT_EQ(kkt_res.lu_passive.size(), dim_passive);
  EXPECT_EQ(kkt_res.lf().size(), dimf);

  kkt_res.Fx.setRandom();
  kkt_res.lx.setRandom();
  EXPECT_TRUE(kkt_res.Fx.head(dimv).isApprox(kkt_res.Fq()));
  EXPECT_TRUE(kkt_res.Fx.tail(dimv).isApprox(kkt_res.Fv()));
  EXPECT_TRUE(kkt_res.lx.head(dimv).isApprox(kkt_res.lq()));
  EXPECT_TRUE(kkt_res.lx.tail(dimv).isApprox(kkt_res.lv()));
}


void SplitKKTResidualTest::testIsApprox(const Robot& robot, 
                                        const ContactStatus& contact_status,
                                        const ImpulseStatus& impulse_status) {
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_res.Fx.setRandom();
  kkt_res.P().setRandom();
  kkt_res.lx.setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu.setRandom();
  kkt_res.lf().setRandom();
  kkt_res.lu_passive.setRandom();

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
  kkt_res.la.setRandom();
  EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.la = kkt_res_ref.la;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.lu.setRandom();
  EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
  kkt_res.lu = kkt_res_ref.lu;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  if (contact_status.hasActiveContacts()) {
    kkt_res_ref.lf().setRandom();
    EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
    kkt_res.lf() = kkt_res_ref.lf();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  else {
    kkt_res_ref.lf().setRandom();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  if (robot.hasFloatingBase()) {
    kkt_res_ref.lu_passive.setRandom();
    EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
    kkt_res.lu_passive = kkt_res_ref.lu_passive;
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  else {
    kkt_res_ref.lu_passive.setRandom();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  if (impulse_status.hasActiveImpulse()) {
    kkt_res_ref.P().setRandom();
    EXPECT_FALSE(kkt_res.isApprox(kkt_res_ref));
    kkt_res.P() = kkt_res_ref.P();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
  else {
    kkt_res_ref.P().setRandom();
    EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  }
}


TEST_F(SplitKKTResidualTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  contact_status.deactivateContact(0);
  impulse_status.deactivateImpulse(0);
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  impulse_status.deactivateImpulse(0);
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


TEST_F(SplitKKTResidualTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  // Both contact and impulse are inactive
  contact_status.deactivateContacts();
  impulse_status.deactivateImpulses();
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  // Contacts are active and impulse are inactive
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  impulse_status.deactivateImpulses();
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  // Contacts are inactive and impulse are active
  contact_status.deactivateContacts();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  // Both contact and impulse are active
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