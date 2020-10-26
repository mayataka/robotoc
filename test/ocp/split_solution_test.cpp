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
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(SplitSolutionTest, fixed_base) {
  Robot robot(fixed_base_urdf_);
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
  EXPECT_TRUE(s.u_passive.size() == 0);
  EXPECT_TRUE(s.nu_passive.size() == 0);
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
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == 0);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 0);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.u_passive.size() == 0);
  EXPECT_TRUE(s_random.nu_passive.size() == 0);
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
}


TEST_F(SplitSolutionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 1);
  EXPECT_TRUE(s.mu[0].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 1);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 0);
  EXPECT_TRUE(s.nu_passive.size() == 0);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 1);
  EXPECT_TRUE(s.mu[0].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 1);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 0);
  EXPECT_TRUE(s.nu_passive.size() == 0);
  EXPECT_TRUE(s.mu_stack().size() == 3);
  EXPECT_TRUE(s.f_stack().size() == 3);
  EXPECT_EQ(s.dimf(), 3);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu0 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd f0 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu[0] = mu0;
  s.a = a;
  s.f[0] = f0;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.set_mu_stack();
  s.set_f_stack();
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.f[0].isApprox(f0));
  EXPECT_TRUE(s.f_stack().isApprox(f0));
  EXPECT_TRUE(s.mu[0].isApprox(mu0));
  EXPECT_TRUE(s.mu_stack().isApprox(mu0));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  s.mu_stack().setZero();
  s.f_stack().setZero();
  s.set_mu();
  s.set_f();
  EXPECT_TRUE(s.mu[0].isZero());
  EXPECT_TRUE(s.f[0].isZero());
  SplitSolution s_random = SplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == 1);
  EXPECT_TRUE(s_random.mu[0].size() == 3);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 1);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.u_passive.size() == 0);
  EXPECT_TRUE(s_random.nu_passive.size() == 0);
  EXPECT_TRUE(s_random.mu_stack().size() == 3);
  EXPECT_TRUE(s_random.f_stack().size() == 3);
  EXPECT_EQ(s_random.dimf(), 3);
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.mu[0].isZero());
  EXPECT_TRUE(s_random.mu[0].isApprox(s_random.mu_stack()));
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.f[0].isZero());
  EXPECT_TRUE(s_random.f[0].isApprox(s_random.f_stack()));
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
}


TEST_F(SplitSolutionTest, floating_base) {
  Robot robot(floating_base_urdf_);
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
  SplitSolution s_random = SplitSolution::Random(robot);
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
  EXPECT_FALSE(s_random.u_passive.isZero());
  EXPECT_FALSE(s_random.nu_passive.isZero());
}


TEST_F(SplitSolutionTest, floating_base_contacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 4);
  EXPECT_TRUE(s.mu[0].size() == 3);
  EXPECT_TRUE(s.mu[1].size() == 3);
  EXPECT_TRUE(s.mu[2].size() == 3);
  EXPECT_TRUE(s.mu[3].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 4);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.f[1].size() == 3);
  EXPECT_TRUE(s.f[2].size() == 3);
  EXPECT_TRUE(s.f[3].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  is_contact_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu.size() == 4);
  EXPECT_TRUE(s.mu[0].size() == 3);
  EXPECT_TRUE(s.mu[1].size() == 3);
  EXPECT_TRUE(s.mu[2].size() == 3);
  EXPECT_TRUE(s.mu[3].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 4);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.f[1].size() == 3);
  EXPECT_TRUE(s.f[2].size() == 3);
  EXPECT_TRUE(s.f[3].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimu());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.u_passive.size() == 6);
  EXPECT_TRUE(s.nu_passive.size() == 6);
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu0 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd mu1 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd mu2 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd mu3 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd f0 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd f1 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd f2 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd f3 = Eigen::VectorXd::Random(3);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu[0] = mu0;
  s.mu[1] = mu1;
  s.mu[2] = mu2;
  s.mu[3] = mu3;
  s.a = a;
  s.f[0] = f0;
  s.f[1] = f1;
  s.f[2] = f2;
  s.f[3] = f3;
  s.q = q;
  s.v = v;
  s.u = u;
  s.beta = beta;
  s.set_mu_stack();
  s.set_f_stack();
  Eigen::VectorXd f_stack_ref = Eigen::VectorXd::Zero(contact_status.dimf());
  Eigen::VectorXd mu_stack_ref = Eigen::VectorXd::Zero(contact_status.dimf());
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      f_stack_ref.segment<3>(dimf_stack) = s.f[i];
      mu_stack_ref.segment<3>(dimf_stack) = s.mu[i];
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.f[0].isApprox(f0));
  EXPECT_TRUE(s.f[1].isApprox(f1));
  EXPECT_TRUE(s.f[2].isApprox(f2));
  EXPECT_TRUE(s.f[3].isApprox(f3));
  EXPECT_TRUE(s.f_stack().isApprox(f_stack_ref));
  EXPECT_TRUE(s.mu[0].isApprox(mu0));
  EXPECT_TRUE(s.mu[1].isApprox(mu1));
  EXPECT_TRUE(s.mu[2].isApprox(mu2));
  EXPECT_TRUE(s.mu[3].isApprox(mu3));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack_ref));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.u.isApprox(u));
  EXPECT_TRUE(s.beta.isApprox(beta));
  s.mu_stack().setZero();
  s.f_stack().setZero();
  s.set_mu();
  s.set_f();
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      EXPECT_TRUE(s.mu[i].isZero());
      EXPECT_TRUE(s.f[i].isZero());
    }
  }
  SplitSolution s_random = SplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu.size() == 4);
  EXPECT_TRUE(s_random.mu[0].size() == 3);
  EXPECT_TRUE(s_random.mu[1].size() == 3);
  EXPECT_TRUE(s_random.mu[2].size() == 3);
  EXPECT_TRUE(s_random.mu[3].size() == 3);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 4);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f[1].size() == 3);
  EXPECT_TRUE(s_random.f[2].size() == 3);
  EXPECT_TRUE(s_random.f[3].size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimu());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.u_passive.size() == 6);
  EXPECT_TRUE(s_random.nu_passive.size() == 6);
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
  EXPECT_EQ(s_random.dimf(), contact_status.dimf());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      EXPECT_FALSE(s_random.mu[i].isZero());
    }
  }
  EXPECT_FALSE(s_random.a.isZero());
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      EXPECT_FALSE(s_random.f[i].isZero());
    }
  }
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}