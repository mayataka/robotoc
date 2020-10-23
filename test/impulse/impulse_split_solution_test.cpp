#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"


namespace idocp {

class ImpulseSplitSolutionTest : public ::testing::Test {
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


TEST_F(ImpulseSplitSolutionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact_velocity.size() == 1);
  EXPECT_TRUE(s.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 1);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimv());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact_velocity.size() == 1);
  EXPECT_TRUE(s.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 1);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == 3);
  EXPECT_TRUE(s.f_stack().size() == 3);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(s.dimc());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f = Eigen::Vector3d::Random(s.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu_contact_velocity[0] = mu_stack;
  s.dv = dv;
  s.f[0] = f;
  s.q = q;
  s.v = v;
  s.beta = beta;
  s.set_mu_stack();
  s.set_f_stack();
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.dv.isApprox(dv));
  EXPECT_TRUE(s.f[0].isApprox(f));
  EXPECT_TRUE(s.f_stack().isApprox(f));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.beta.isApprox(beta));
  s.mu_stack().setZero();
  s.f_stack().setZero();
  s.set_mu_contact();
  s.set_f();
  EXPECT_TRUE(s.f[0].isZero());
  EXPECT_EQ(s.dimf(), 3);
  EXPECT_EQ(s.dimc(), 3);
  ImpulseSplitSolution s_random = ImpulseSplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact_velocity.size() == 1);
  EXPECT_TRUE(s_random.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s_random.mu_stack().size() == 3);
  EXPECT_TRUE(s_random.dv.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 1);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f_stack().size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
}


TEST_F(ImpulseSplitSolutionTest, floating_base_contacts_full_contacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true, true, true, true};
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact_velocity.size() == 4);
  EXPECT_TRUE(s.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[1].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[2].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[3].size() == 3);
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 4);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.f[1].size() == 3);
  EXPECT_TRUE(s.f[2].size() == 3);
  EXPECT_TRUE(s.f[3].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(s.dimc());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f = Eigen::VectorXd::Random(s.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.mu_contact_velocity[0] = mu_stack.segment<3>(0);
  s.mu_contact_velocity[1] = mu_stack.segment<3>(3);
  s.mu_contact_velocity[2] = mu_stack.segment<3>(6);
  s.mu_contact_velocity[3] = mu_stack.segment<3>(9);
  s.dv = dv;
  s.f[0] = f.segment<3>(0);
  s.f[1] = f.segment<3>(3);
  s.f[2] = f.segment<3>(6);
  s.f[3] = f.segment<3>(9);
  s.q = q;
  s.v = v;
  s.beta = beta;
  s.set_mu_stack();
  s.set_f_stack();
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.dv.isApprox(dv));
  EXPECT_TRUE(s.f_stack().isApprox(f));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.beta.isApprox(beta));
  s.mu_stack().setZero();
  s.f_stack().setZero();
  s.set_mu_contact();
  s.set_f();
  EXPECT_TRUE(s.mu_contact_velocity[0].isZero());
  EXPECT_TRUE(s.mu_contact_velocity[1].isZero());
  EXPECT_TRUE(s.mu_contact_velocity[2].isZero());
  EXPECT_TRUE(s.mu_contact_velocity[3].isZero());
  EXPECT_TRUE(s.f[0].isZero());
  EXPECT_TRUE(s.f[1].isZero());
  EXPECT_TRUE(s.f[2].isZero());
  EXPECT_TRUE(s.f[3].isZero());
  EXPECT_EQ(s.dimf(), contact_status.dimf());
  EXPECT_EQ(s.dimc(), contact_status.dimf());
  ImpulseSplitSolution s_random = ImpulseSplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact_velocity.size() == 4);
  EXPECT_TRUE(s_random.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[1].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[2].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[3].size() == 3);
  EXPECT_TRUE(s_random.dv.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 4);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f[1].size() == 3);
  EXPECT_TRUE(s_random.f[2].size() == 3);
  EXPECT_TRUE(s_random.f[3].size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
}


TEST_F(ImpulseSplitSolutionTest, floating_base_contacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  ContactStatus contact_status(is_contact_active.size());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitSolution s(robot);
  s.setContactStatus(contact_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.mu_contact_velocity.size() == 4);
  EXPECT_TRUE(s.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[1].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[2].size() == 3);
  EXPECT_TRUE(s.mu_contact_velocity[3].size() == 3);
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 4);
  EXPECT_TRUE(s.f[0].size() == 3);
  EXPECT_TRUE(s.f[1].size() == 3);
  EXPECT_TRUE(s.f[2].size() == 3);
  EXPECT_TRUE(s.f[3].size() == 3);
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s.f_stack().size() == contact_status.dimf());
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(s.dimc());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f = Eigen::VectorXd::Random(s.dimf());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  ImpulseSplitSolution s_random = ImpulseSplitSolution::Random(robot, contact_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_contact_velocity.size() == 4);
  EXPECT_TRUE(s_random.mu_contact_velocity[0].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[1].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[2].size() == 3);
  EXPECT_TRUE(s_random.mu_contact_velocity[3].size() == 3);
  EXPECT_TRUE(s_random.dv.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 4);
  EXPECT_TRUE(s_random.f[0].size() == 3);
  EXPECT_TRUE(s_random.f[1].size() == 3);
  EXPECT_TRUE(s_random.f[2].size() == 3);
  EXPECT_TRUE(s_random.f[3].size() == 3);
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.mu_stack().size() == contact_status.dimf());
  EXPECT_TRUE(s_random.f_stack().size() == contact_status.dimf());
}



} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}