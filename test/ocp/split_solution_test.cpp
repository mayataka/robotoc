#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
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
  EXPECT_TRUE(s.mu_contact.size() == 0);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.q.size() == robot.dimv());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.mu_floating_base().size() == 0);
  EXPECT_TRUE(s.mu_contacts().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  s.setContactStatus(robot);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.a = a;
  s.q = q;
  s.v = v;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_EQ(s.dimc(), 0);
  SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact.size() == 0);
  EXPECT_TRUE(s_random.mu_stack().size() == 0);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 0);
  EXPECT_TRUE(s_random.f_stack().size() == 0);
  EXPECT_TRUE(s_random.q.size() == robot.dimv());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_EQ(s_random.dimc(), 0);
}


TEST_F(SplitSolutionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  Robot robot(fixed_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status = {false};
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact.size() == 1);
  EXPECT_TRUE(s.mu_contact[0].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 1);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimv());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.mu_floating_base().size() == 0);
  EXPECT_TRUE(s.mu_contacts().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  contact_status = {true};
  robot.setContactStatus(contact_status);
  s.setContactStatus(robot);
  EXPECT_TRUE(s.mu_stack().size() == 3);
  EXPECT_TRUE(s.mu_floating_base().size() == 0);
  EXPECT_TRUE(s.mu_contacts().size() == 3);
  EXPECT_TRUE(s.f_stack().size() == 3);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(robot.dimf()+robot.dim_passive());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd f = Eigen::Vector3d::Random();
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu_contact[0] = mu_stack;
  s.a = a;
  s.f[0] = f;
  s.q = q;
  s.v = v;
  s.set_mu_stack();
  s.set_f_stack();
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu_contacts().isApprox(mu_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.f[0].isApprox(f));
  EXPECT_TRUE(s.f_stack().isApprox(f));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  s.mu_stack().setZero();
  s.f_stack().setZero();
  s.set_mu_contact();
  s.set_f();
  EXPECT_TRUE(s.mu_contacts().isZero());
  EXPECT_TRUE(s.f[0].isZero());
  EXPECT_EQ(s.dimf(), 3);
  EXPECT_EQ(s.dimc(), 3);
  SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact.size() == 1);
  EXPECT_TRUE(s_random.mu_contact[0].size() == 3);
  EXPECT_TRUE(s_random.mu_stack().size() == 3);
  EXPECT_TRUE(s_random.mu_floating_base().size() == 0);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 1);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f_stack().size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimv());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.mu_stack().isZero());
  EXPECT_FALSE(s_random.mu_contact[0].isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.f[0].isZero());
  EXPECT_FALSE(s_random.f_stack().isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_EQ(s_random.dimf(), 3);
  EXPECT_EQ(s_random.dimc(), 3);
}


TEST_F(SplitSolutionTest, floating_base) {
  Robot robot(floating_base_urdf_);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact.size() == 0);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 6);
  EXPECT_TRUE(s.mu_floating_base().size() == 6);
  EXPECT_TRUE(s.mu_contacts().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(robot.dimf()+robot.dim_passive());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu_floating_base() = mu_stack;
  s.a = a;
  s.q = q;
  s.v = v;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu_floating_base().isApprox(mu_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  s.mu_stack().setZero();
  EXPECT_TRUE(s.mu_floating_base().isZero());
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_EQ(s.dimc(), 6);
  SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact.size() == 0);
  EXPECT_TRUE(s_random.mu_stack().size() == 6);
  EXPECT_TRUE(s_random.mu_floating_base().size() == 6);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 0);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.mu_stack().isZero());
  EXPECT_FALSE(s_random.mu_floating_base().isZero());
  EXPECT_FALSE(s_random.a.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_EQ(s_random.dimc(), 6);
}


TEST_F(SplitSolutionTest, floating_base_contacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  Robot robot(floating_base_urdf_, contact_frames, mu, 0, 0);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitSolution s(robot);
  s.setContactStatus(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact.size() == 4);
  EXPECT_TRUE(s.mu_contact[0].size() == 3);
  EXPECT_TRUE(s.mu_contact[1].size() == 3);
  EXPECT_TRUE(s.mu_contact[2].size() == 3);
  EXPECT_TRUE(s.mu_contact[3].size() == 3);
  EXPECT_TRUE(s.a.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 4);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.f[1].size() == 3);
  EXPECT_TRUE(s.f[2].size() == 3);
  EXPECT_TRUE(s.f[3].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.u.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 6+robot.dimf());
  EXPECT_TRUE(s.mu_floating_base().size() == 6);
  EXPECT_TRUE(s.mu_contacts().size() == robot.dimf());
  EXPECT_TRUE(s.f_stack().size() == robot.dimf());
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(robot.dimf()+robot.dim_passive());
  const Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd f_stack = Eigen::VectorXd::Random(robot.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu_stack() = mu_stack;
  s.a = a;
  s.f_stack() = f_stack;
  s.q = q;
  s.v = v;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.mu_floating_base().isApprox(mu_stack.head(6)));
  EXPECT_TRUE(s.mu_contacts().isApprox(mu_stack.tail(robot.dimf())));
  EXPECT_TRUE(s.a.isApprox(a));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  s.set_mu_contact();
  s.set_f();
  int sum = 0;
  for (int i=0; i<4; ++i) {
    if (robot.is_contact_active(i)) {
      EXPECT_TRUE(s.mu_contact[i].isApprox(mu_stack.template segment<3>(6+3*sum)));
      ++sum;
    }
  }
  sum = 0;
  for (int i=0; i<4; ++i) {
    if (robot.is_contact_active(i)) {
      EXPECT_TRUE(s.f[i].isApprox(f_stack.template segment<3>(3*sum)));
      ++sum;
    }
  }
  std::vector<Eigen::Vector3d> f_ref, mu_ref;
  for (int i=0; i<4; ++i) {
    f_ref.push_back(Eigen::Vector3d::Random());
    mu_ref.push_back(Eigen::Vector3d::Random());
  }
  for (int i=0; i<4; ++i) {
    s.f[i] = f_ref[i];
    s.mu_contact[i] = mu_ref[i];
  }
  s.set_f_stack();
  s.set_mu_stack();
  sum = 0;
  EXPECT_TRUE(mu_stack.head(6).isApprox(s.mu_stack().head(6)));
  EXPECT_TRUE(mu_stack.head(6).isApprox(s.mu_floating_base()));
  for (int i=0; i<4; ++i) {
    if (robot.is_contact_active(i)) {
      EXPECT_TRUE(mu_ref[i].isApprox(s.mu_stack().template segment<3>(6+3*sum)));
      ++sum;
    }
  }
  sum = 0;
  for (int i=0; i<4; ++i) {
    if (robot.is_contact_active(i)) {
      EXPECT_TRUE(f_ref[i].isApprox(s.f_stack().template segment<3>(3*sum)));
      ++sum;
    }
  }
  EXPECT_EQ(s.dimf(), robot.dimf());
  EXPECT_EQ(s.dimc(), 6+robot.dimf());
  SplitSolution s_random = SplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact.size() == 4);
  EXPECT_TRUE(s_random.mu_contact[0].size() == 3);
  EXPECT_TRUE(s_random.mu_contact[1].size() == 3);
  EXPECT_TRUE(s_random.mu_contact[2].size() == 3);
  EXPECT_TRUE(s_random.mu_contact[3].size() == 3);
  EXPECT_TRUE(s_random.mu_stack().size() == s.dimc());
  EXPECT_TRUE(s_random.mu_floating_base().size() == 6);
  EXPECT_TRUE(s_random.a.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 4);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f[1].size() == 3);
  EXPECT_TRUE(s_random.f[2].size() == 3);
  EXPECT_TRUE(s_random.f[3].size() == 3);
  EXPECT_TRUE(s_random.f_stack().size() == s.dimf());
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.u.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.mu_stack().isZero());
  EXPECT_FALSE(s_random.mu_floating_base().isZero());
  EXPECT_FALSE(s_random.a.isZero());
  if (s_random.dimf() > 0) {
    EXPECT_FALSE(s_random.f_stack().isZero());
  }
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.u.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  EXPECT_EQ(s_random.dimf(), s.dimf());
  EXPECT_EQ(s_random.dimc(), s.dimc());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}