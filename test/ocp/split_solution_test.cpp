#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"


namespace idocp {

class SplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void TestWithoutContacts(const Robot& robot);
  static void TestWithContacts(const Robot& robot, const ContactStatus& contact_status);
  static void TestIsApprox(const Robot& robot, const ContactStatus& contact_status);
  static void TestIntegrate(const Robot& robot, const ContactStatus& contact_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitSolutionTest::TestWithoutContacts(const Robot& robot) { 
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
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_passive = Eigen::VectorXd::Random(6);
  const Eigen::VectorXd nu_passive = Eigen::VectorXd::Random(6);
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.u_passive = u_passive;
  s.nu_passive = nu_passive;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.u_passive.isApprox(u_passive));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  s.u_passive.setZero();
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
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (robot.has_floating_base()) {
    EXPECT_FALSE(s.u_passive.isZero());
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.u_passive.isZero());
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
  EXPECT_TRUE(s_random.u_passive.size() == 6);
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == 0);
  EXPECT_TRUE(s_random.f_stack().size() == 0);
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  if (robot.has_floating_base()) {
    EXPECT_FALSE(s_random.u_passive.isZero());
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.u_passive.isZero());
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::TestWithContacts(const Robot& robot, 
                                         const ContactStatus& contact_status) { 
  std::random_device rnd;
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_EQ(s.isContactActive(i), contact_status.isContactActive(i));
  }
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u_passive = Eigen::VectorXd::Random(6);
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
  s.u_passive = u_passive;
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
  EXPECT_TRUE(s.u_passive.isApprox(u_passive));
  EXPECT_TRUE(s.nu_passive.isApprox(nu_passive));
  EXPECT_TRUE(s.f_stack().isApprox(f_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  s.set_f_vector();
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_mu_vector();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
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
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.u_passive.setZero();
  s.nu_passive.setZero();
  s = SplitSolution(robot);
  s.setRandom(robot, contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
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
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isContactActive(i)) {
      EXPECT_FALSE(s.mu[i].isZero());
    }
  }
  if (robot.has_floating_base()) {
    EXPECT_FALSE(s.u_passive.isZero());
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s.u_passive.isZero());
    EXPECT_TRUE(s.nu_passive.isZero());
  }
  const SplitSolution s_random = SplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s_random.mu[i].size() == 3);
  }
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s_random.f[i].size() == 3);
  }
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.u_passive.size() == 6);
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s_random.dimf(), contact_status.dimf());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
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
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s_random.isContactActive(i)) {
      EXPECT_FALSE(s_random.mu[i].isZero());
    }
  }
  if (robot.has_floating_base()) {
    EXPECT_FALSE(s_random.u_passive.isZero());
    EXPECT_FALSE(s_random.nu_passive.isZero());
  }
  else {
    EXPECT_TRUE(s_random.u_passive.isZero());
    EXPECT_TRUE(s_random.nu_passive.isZero());
  }
}


void SplitSolutionTest::TestIsApprox(const Robot& robot, 
                                     const ContactStatus& contact_status) {
  SplitSolution s(robot);
  s.setRandom(robot, contact_status);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.a.isZero());
  EXPECT_FALSE(s.u.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (robot.has_floating_base()) {
    EXPECT_FALSE(s.u_passive.isZero());
    EXPECT_FALSE(s.nu_passive.isZero());
  }
  if (s.hasActiveContacts()) {
    EXPECT_FALSE(s.f_stack().isZero());
    EXPECT_FALSE(s.mu_stack().isZero());
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
  if (robot.has_floating_base()) {
    s_ref.u_passive.setRandom();
    s_ref.nu_passive.setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.u_passive = s.u_passive;
    s_ref.nu_passive = s.nu_passive;
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  else {
    s_ref.u_passive.setRandom();
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
}


void SplitSolutionTest::TestIntegrate(const Robot& robot, 
                                      const ContactStatus& contact_status) {
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
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
  if (contact_status.hasActiveContacts()) {
    s_ref.f_stack().noalias() += step_size * d.df();
    s_ref.set_f_vector();
    s_ref.mu_stack().noalias() += step_size * d.dmu();
    s_ref.set_mu_vector();
  }
  if (robot.has_floating_base()) {
    s_ref.u_passive.noalias() += step_size * d.du_passive;
    s_ref.nu_passive.noalias() += step_size * d.dnu_passive;
  }
  EXPECT_TRUE(s.isApprox(s_ref));
}


TEST_F(SplitSolutionTest, fixedBase) {
  Robot robot_without_contacts(fixed_base_urdf);
  TestWithoutContacts(robot_without_contacts);
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  TestIsApprox(robot, contact_status);
  TestIntegrate(robot, contact_status);
  contact_status.activateContact(0);
  TestWithContacts(robot, contact_status);
  TestIsApprox(robot, contact_status);
  TestIntegrate(robot, contact_status);
}


TEST_F(SplitSolutionTest, floatingBase) {
  Robot robot_without_contacts(floating_base_urdf);
  TestWithoutContacts(robot_without_contacts);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  TestIsApprox(robot, contact_status);
  TestIntegrate(robot, contact_status);
  std::random_device rnd;
  is_contact_active.clear();
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  TestWithContacts(robot, contact_status);
  TestIsApprox(robot, contact_status);
  TestIntegrate(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}