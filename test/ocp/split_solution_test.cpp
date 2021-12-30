#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/ocp/split_solution.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class SplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot);
  static void test(const Robot& robot, const ContactStatus& contact_status);
  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void test(const Robot& robot, const ContactStatus& contact_status, 
                   const ImpulseStatus& impulse_status);
  static void test_isApprox(const Robot& robot, 
                           const ContactStatus& contact_status, 
                           const ImpulseStatus& impulse_status);
  static void test_integrate(const Robot& robot, 
                            const ContactStatus& contact_status, 
                            const ImpulseStatus& impulse_status);

  double dt;
};


void SplitSolutionTest::test(const Robot& robot) { 
  SplitSolution s(robot);
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.u.size(), robot.dimu());
  EXPECT_EQ(s.a.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxNumContacts());
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.mu.size(), robot.maxNumContacts());
  EXPECT_EQ(s.nu_passive.size(), robot.dim_passive());
  EXPECT_EQ(s.f_stack().size(), 0);
  EXPECT_EQ(s.mu_stack().size(), 0);
  EXPECT_EQ(s.xi_stack().size(), 0);
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_EQ(s.dimi(), 0);
}


void SplitSolutionTest::test(const Robot& robot, 
                             const ContactStatus& contact_status) { 
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.u.size(), robot.dimu());
  EXPECT_EQ(s.a.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxNumContacts());
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.mu.size(), robot.maxNumContacts());
  EXPECT_EQ(s.nu_passive.size(), robot.dim_passive());
  EXPECT_EQ(s.f_stack().size(), contact_status.dimf());
  EXPECT_EQ(s.mu_stack().size(), contact_status.dimf());
  EXPECT_EQ(s.xi_stack().size(), 0);
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  EXPECT_EQ(s.dimi(), 0);
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
    EXPECT_EQ(s.isContactActive()[i], contact_status.isContactActive(i));
  }

  EXPECT_NO_THROW(
    std::cout << s << std::endl;
  );
}


void SplitSolutionTest::test(const Robot& robot, 
                             const ImpulseStatus& impulse_status) {
  SplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.u.size(), robot.dimu());
  EXPECT_EQ(s.a.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxNumContacts());
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.mu.size(), robot.maxNumContacts());
  EXPECT_EQ(s.nu_passive.size(), robot.dim_passive());
  EXPECT_EQ(s.f_stack().size(), 0);
  EXPECT_EQ(s.mu_stack().size(), 0);
  EXPECT_EQ(s.xi_stack().size(), impulse_status.dimi());
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_EQ(s.dimi(), impulse_status.dimi());

  EXPECT_NO_THROW(
    std::cout << s << std::endl;
  );
}


void SplitSolutionTest::test(const Robot& robot, 
                             const ContactStatus& contact_status,
                             const ImpulseStatus& impulse_status) { 
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  s.setImpulseStatus(impulse_status);
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.u.size(), robot.dimu());
  EXPECT_EQ(s.a.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxNumContacts());
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.mu.size(), robot.maxNumContacts());
  EXPECT_EQ(s.nu_passive.size(), robot.dim_passive());
  EXPECT_EQ(s.f_stack().size(), contact_status.dimf());
  EXPECT_EQ(s.mu_stack().size(), contact_status.dimf());
  EXPECT_EQ(s.xi_stack().size(), impulse_status.dimi());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimi());
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
    EXPECT_EQ(s.isContactActive()[i], contact_status.isContactActive(i));
  }

  EXPECT_NO_THROW(
    std::cout << s << std::endl;
  );
}


void SplitSolutionTest::test_isApprox(const Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const ImpulseStatus& impulse_status) {
  auto s = SplitSolution::Random(robot, contact_status, impulse_status);
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (s.hasActiveContacts()) {
    EXPECT_FALSE(s.f_stack().isZero());
    EXPECT_FALSE(s.mu_stack().isZero());
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  if (s.hasActiveImpulse()) {
    EXPECT_FALSE(s.xi_stack().isZero());
  }
  auto s_ref = s;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.q.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.q = s.q;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.v.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.v = s.v;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.a.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.a = s.a;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.u.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.u = s.u;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.lmd.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.lmd = s.lmd;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.gmm.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.gmm = s.gmm;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.beta.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.beta = s.beta;
  EXPECT_TRUE(s.isApprox(s_ref));
  if (robot.hasFloatingBase()) {
    s_ref.nu_passive.setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.nu_passive = s.nu_passive;
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  else {
    s_ref.nu_passive.setRandom();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  if (s.hasActiveContacts()) {
    s_ref.f_stack().setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.f_stack() = s.f_stack();
    s_ref.set_f_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
    s_ref.mu_stack().setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.mu_stack() = s.mu_stack();
    s_ref.set_mu_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  else {
    s_ref.f_stack().setRandom();
    s_ref.set_f_vector();
    s_ref.mu_stack().setRandom();
    s_ref.set_mu_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  if (s.hasActiveImpulse()) {
    s_ref.xi_stack().setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.xi_stack() = s.xi_stack();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  else {
    s_ref.xi_stack().setRandom();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  s_ref.setRandom(robot);
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.copyPrimal(s);
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.setRandom(robot);
  s_ref.copyDual(s);
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.setRandom(robot);
  s_ref.copyPrimal(s);
  s_ref.copyDual(s);
  EXPECT_TRUE(s.isApprox(s_ref));
}


void SplitSolutionTest::test_integrate(const Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const ImpulseStatus& impulse_status) {
  SplitSolution s = SplitSolution::Random(robot, contact_status, impulse_status);
  const SplitDirection d = SplitDirection::Random(robot, contact_status, impulse_status);
  SplitSolution s_ref = s;
  const double step_size = 0.3;
  s.integrate(robot, step_size, d);
  s_ref.lmd.noalias() += step_size * d.dlmd();
  s_ref.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s_ref.q);
  s_ref.v.noalias() += step_size * d.dv();
  s_ref.a.noalias() += step_size * d.da();
  s_ref.u.noalias() += step_size * d.du;
  s_ref.beta.noalias() += step_size * d.dbeta();
  if (robot.hasFloatingBase()) {
    s_ref.nu_passive.noalias() += step_size * d.dnu_passive;
  }
  if (contact_status.hasActiveContacts()) {
    s_ref.f_stack().noalias() += step_size * d.df();
    s_ref.set_f_vector();
    s_ref.mu_stack().noalias() += step_size * d.dmu();
    s_ref.set_mu_vector();
  }
  if (impulse_status.hasActiveImpulse()) {
    s_ref.xi_stack().noalias() += step_size * d.dxi();
  }
  EXPECT_TRUE(s.isApprox(s_ref));
}


TEST_F(SplitSolutionTest, fixedBase) {
  auto robot_without_contacts = testhelper::CreateRobotManipulator();
  test(robot_without_contacts);
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  test(robot, contact_status);
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
  test(robot, contact_status, impulse_status);
}


TEST_F(SplitSolutionTest, floatingBase) {
  auto robot_without_contacts = testhelper::CreateQuadrupedalRobot();
  test(robot_without_contacts);
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status);
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
  test(robot, contact_status, impulse_status);
}


TEST_F(SplitSolutionTest, humanoidRobot) {
  auto robot_without_contacts = testhelper::CreateHumanoidRobot();
  test(robot_without_contacts);
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test(robot, contact_status);
  test_isApprox(robot, contact_status, impulse_status);
  test_integrate(robot, contact_status, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
  test(robot, contact_status, impulse_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}