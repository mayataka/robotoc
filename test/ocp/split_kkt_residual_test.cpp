#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SplitKKTResidualTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  static void test(const Robot& robot, const ContactStatus& contact_status);
  static void testIsApprox(const Robot& robot, const ContactStatus& contact_status);

  virtual void TearDown() {
  }

  double dt;
};


void SplitKKTResidualTest::test(const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimx = 2 * robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  EXPECT_EQ(kkt_res.dimf(), dimf);
  EXPECT_EQ(kkt_res.Fx.size(), dimx);
  EXPECT_EQ(kkt_res.Fq().size(), dimv);
  EXPECT_EQ(kkt_res.Fv().size(), dimv);
  EXPECT_EQ(kkt_res.lx.size(), dimx);
  EXPECT_EQ(kkt_res.lq().size(), dimv);
  EXPECT_EQ(kkt_res.lv().size(), dimv);
  EXPECT_EQ(kkt_res.la.size(), dimv);
  EXPECT_EQ(kkt_res.lu.size(), dimu);
  EXPECT_EQ(kkt_res.lf().size(), dimf);

  kkt_res.Fx.setRandom();
  kkt_res.lx.setRandom();
  EXPECT_TRUE(kkt_res.Fx.head(dimv).isApprox(kkt_res.Fq()));
  EXPECT_TRUE(kkt_res.Fx.tail(dimv).isApprox(kkt_res.Fv()));
  EXPECT_TRUE(kkt_res.lx.head(dimv).isApprox(kkt_res.lq()));
  EXPECT_TRUE(kkt_res.lx.tail(dimv).isApprox(kkt_res.lv()));

  kkt_res.lu.setRandom();
  kkt_res.la.setRandom();
  kkt_res.lf().setRandom();
  const double nrm = kkt_res.squaredNormKKTResidual();
  const double nrm_ref = kkt_res.Fx.squaredNorm() 
                          + kkt_res.lx.squaredNorm()
                          + kkt_res.lu.squaredNorm()
                          + kkt_res.la.squaredNorm()
                          + kkt_res.lf().squaredNorm();
  EXPECT_DOUBLE_EQ(nrm, nrm_ref);

  const double vio = kkt_res.l1NormConstraintViolation();
  const double vio_ref = kkt_res.Fx.template lpNorm<1>();
  EXPECT_DOUBLE_EQ(vio, vio_ref);
}


void SplitKKTResidualTest::testIsApprox(const Robot& robot, const ContactStatus& contact_status) {
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.Fx.setRandom();
  kkt_res.lx.setRandom();
  kkt_res.la.setRandom();
  kkt_res.lu.setRandom();
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
}


TEST_F(SplitKKTResidualTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  contact_status.deactivateContact(0);
  test(robot, contact_status);
  testIsApprox(robot, contact_status);
  contact_status.activateContact(0);
  test(robot, contact_status);
  testIsApprox(robot, contact_status);
}


TEST_F(SplitKKTResidualTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  contact_status.deactivateContacts();
  test(robot, contact_status);
  testIsApprox(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status);
  testIsApprox(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}