#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_solution.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static void Test(const Robot& robot);
  static void Test(const Robot& robot, const ContactStatus& contact_status);
  static void Test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void Test(const Robot& robot, const ContactStatus& contact_status, 
                   const ImpulseStatus& impulse_status);
  static void TestIsApprox(const Robot& robot, 
                           const ContactStatus& contact_status, 
                           const ImpulseStatus& impulse_status);
  static void TestCopy(const Robot& robot, const ContactStatus& contact_status,
                       const ImpulseStatus& impulse_status);
  static void TestIntegrate(const Robot& robot, 
                            const ContactStatus& contact_status, 
                            const ImpulseStatus& impulse_status);

  double dt;
};


void SplitSolutionTest::Test(const Robot& robot) { 
  std::random_device rnd;
  SplitSolution s(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 0);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd nu_passive = Eigen::VectorXd::Random(6);
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.nu_passive = nu_passive;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  s.nu_passive.setZero();
  s.setRandom(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 0);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.nu_passive.isZero());
  }
  const SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == 0);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 0);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == 0);
  EXPECT_TRUE(s_random.f_stack().size() == 0);
  EXPECT_TRUE(s_random.xi_stack().size() == 0);
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::Test(const Robot& robot, 
                             const ContactStatus& contact_status) { 
  std::random_device rnd;
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
  }
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd nu_passive = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd f_stack = Eigen::VectorXd::Random(contact_status.dimf());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(contact_status.dimf());
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.nu_passive = nu_passive;
  s.f_stack() = f_stack;
  s.mu_stack() = mu_stack;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  EXPECT_TRUE(s.f_stack().isApprox(f_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  s.set_f_vector();
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_mu_vector();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.f_stack().setZero();
  s.mu_stack().setZero();
  s.set_f_stack();
  s.set_mu_stack();
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.nu_passive.setZero();
  s = SplitSolution(robot);
  s.setRandom(robot, contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
  }
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  EXPECT_FALSE(s.f_stack().isZero());
  EXPECT_FALSE(s.mu_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.nu_passive.isZero());
  }
  const SplitSolution s_random = SplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.mu[i].size() == 3);
  }
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.f[i].size() == 3);
  }
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s_random.dimf(), contact_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s_random.isContactActive(i), contact_status.isContactActive(i));
  }
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_FALSE(s_random.f_stack().isZero());
  EXPECT_FALSE(s_random.mu_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s_random.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::Test(const Robot& robot, 
                             const ImpulseStatus& impulse_status) {
  std::random_device rnd;
  SplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd nu_passive = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd xi_stack = Eigen::VectorXd::Random(impulse_status.dimf());
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.nu_passive = nu_passive;
  s.xi_stack() = xi_stack;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  EXPECT_TRUE(s.xi_stack().isApprox(xi_stack));
  s.xi_stack().setZero();
  s.nu_passive.setZero();
  s = SplitSolution(robot);
  s.setRandom(robot, impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimf());
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  EXPECT_FALSE(s.xi_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.nu_passive.isZero());
  }
  s.setImpulseStatus();
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimi(), 0);
  const SplitSolution s_random = SplitSolution::Random(robot, impulse_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.mu[i].size() == 3);
  }
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.f[i].size() == 3);
  }
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == 0);
  EXPECT_TRUE(s_random.f_stack().size() == 0);
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_EQ(s_random.dimi(), impulse_status.dimf());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_TRUE(s_random.f_stack().isZero());
  EXPECT_TRUE(s_random.mu_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_TRUE(s_random.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::Test(const Robot& robot, 
                             const ContactStatus& contact_status,
                             const ImpulseStatus& impulse_status) { 
  std::random_device rnd;
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  s.setImpulseStatus(impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimf());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
  }
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd nu_passive = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd f_stack = Eigen::VectorXd::Random(contact_status.dimf());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(contact_status.dimf());
  const Eigen::VectorXd xi_stack = Eigen::VectorXd::Random(impulse_status.dimf());
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.nu_passive = nu_passive;
  s.f_stack() = f_stack;
  s.mu_stack() = mu_stack;
  s.xi_stack() = xi_stack;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  EXPECT_TRUE(s.f_stack().isApprox(f_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.xi_stack().isApprox(xi_stack));
  s.set_f_vector();
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_mu_vector();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.f_stack().setZero();
  s.mu_stack().setZero();
  s.set_f_stack();
  s.set_mu_stack();
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.nu_passive.setZero();
  s = SplitSolution(robot);
  s.setRandom(robot, contact_status, impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimf());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
  }
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  EXPECT_FALSE(s.f_stack().isZero());
  EXPECT_FALSE(s.mu_stack().isZero());
  EXPECT_FALSE(s.xi_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.nu_passive.isZero());
  }
  s.setImpulseStatus();
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimi(), 0);
  const SplitSolution s_random = SplitSolution::Random(robot, contact_status, impulse_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.mu[i].size() == 3);
  }
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == robot.maxPointContacts());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s_random.f[i].size() == 3);
  }
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.xi_stack().size() == impulse_status.dimf());
  EXPECT_EQ(s_random.dimf(), contact_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s_random.isContactActive(i), contact_status.isContactActive(i));
  }
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_FALSE(s_random.f_stack().isZero());
  EXPECT_FALSE(s_random.mu_stack().isZero());
  EXPECT_FALSE(s_random.xi_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s_random.mu[i].isZero());
    }
  }
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::TestIsApprox(const Robot& robot, 
                                     const ContactStatus& contact_status,
                                     const ImpulseStatus& impulse_status) {
  SplitSolution s(robot);
  s.setRandom(robot, contact_status, impulse_status);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (robot.hasFloatingBase()) {
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  if (s.hasActiveContacts()) {
    EXPECT_FALSE(s.f_stack().isZero());
    EXPECT_FALSE(s.mu_stack().isZero());
  }
  if (s.hasActiveImpulse()) {
    EXPECT_FALSE(s.xi_stack().isZero());
  }
  SplitSolution s_ref = s;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.lmd.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.lmd = s.lmd;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.gmm.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.gmm = s.gmm;
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
}


void SplitSolutionTest::TestCopy(const Robot& robot, 
                                 const ContactStatus& contact_status,
                                 const ImpulseStatus& impulse_status) {
  SplitSolution s(robot);
  s.setRandom(robot, contact_status, impulse_status);
  SplitSolution s_new(robot);
  s_new.copy(s);
  EXPECT_TRUE(s_new.isApprox(s));
}


void SplitSolutionTest::TestIntegrate(const Robot& robot, 
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
  s_ref.u.noalias() += step_size * d.du();
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
  auto robot_without_contacts = testhelper::CreateFixedBaseRobot();
  Test(robot_without_contacts);
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  TestIsApprox(robot, contact_status, impulse_status);
  TestCopy(robot, contact_status, impulse_status);
  TestIntegrate(robot, contact_status, impulse_status);
  contact_status.activateContact(0);
  Test(robot, contact_status);
  TestIsApprox(robot, contact_status, impulse_status);
  TestCopy(robot, contact_status, impulse_status);
  TestIntegrate(robot, contact_status, impulse_status);
  impulse_status.activateImpulse(0);
  Test(robot, impulse_status);
  Test(robot, contact_status, impulse_status);
}


TEST_F(SplitSolutionTest, floatingBase) {
  auto robot_without_contacts = testhelper::CreateFloatingBaseRobot();
  Test(robot_without_contacts);
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  auto impulse_status = robot.createImpulseStatus();
  TestIsApprox(robot, contact_status, impulse_status);
  TestCopy(robot, contact_status, impulse_status);
  TestIntegrate(robot, contact_status, impulse_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  Test(robot, contact_status);
  TestIsApprox(robot, contact_status, impulse_status);
  TestCopy(robot, contact_status, impulse_status);
  TestIntegrate(robot, contact_status, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  Test(robot, impulse_status);
  Test(robot, contact_status, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}