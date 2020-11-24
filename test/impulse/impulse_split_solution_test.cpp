#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"


namespace idocp {

class ImpulseSplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void TestWithoutImpulses(const Robot& robot);
  static void TestWithImpulses(const Robot& robot, const ImpulseStatus& impulse_status);
  static void TestIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);
  static void TestIntegrate(const Robot& robot, const ImpulseStatus& impulse_status);

  double dtau_;
  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseSplitSolutionTest::TestWithoutImpulses(const Robot& robot) { 
  std::random_device rnd;
  ImpulseSplitSolution s(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.mu.size() == 0);
  EXPECT_TRUE(s.xi.size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  s.lmd = lmd;
  s.gmm = gmm;
  s.dv = dv;
  s.q = q;
  s.v = v;
  s.dv = dv;
  s.beta = beta;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.dv.isApprox(dv));
  EXPECT_TRUE(s.beta.isApprox(beta));
  s.setRandom(robot);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == 0);
  EXPECT_TRUE(s.mu.size() == 0);
  EXPECT_TRUE(s.xi.size() == 0);
  EXPECT_TRUE(s.f_stack().size() == 0);
  EXPECT_TRUE(s.mu_stack().size() == 0);
  EXPECT_TRUE(s.xi_stack().size() == 0);
  EXPECT_EQ(s.dimf(), 0);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.dv.isZero());
  EXPECT_FALSE(s.beta.isZero());
  const ImpulseSplitSolution s_random = ImpulseSplitSolution::Random(robot);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.dv.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == 0);
  EXPECT_TRUE(s_random.mu.size() == 0);
  EXPECT_TRUE(s_random.xi.size() == 0);
  EXPECT_TRUE(s_random.f_stack().size() == 0);
  EXPECT_TRUE(s_random.mu_stack().size() == 0);
  EXPECT_TRUE(s_random.xi_stack().size() == 0);
  EXPECT_EQ(s_random.dimf(), 0);
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.dv.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
}


void ImpulseSplitSolutionTest::TestWithImpulses(const Robot& robot, 
                                                const ImpulseStatus& impulse_status) { 
  std::random_device rnd;
  ImpulseSplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.max_point_contacts());
  EXPECT_TRUE(s.mu.size() == robot.max_point_contacts());
  EXPECT_TRUE(s.xi.size() == robot.max_point_contacts());
  EXPECT_TRUE(s.f_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s.mu_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimp());
  EXPECT_EQ(s.dimf(), impulse_status.dimp());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
    EXPECT_TRUE(s.mu[i].size() == 3);
    EXPECT_TRUE(s.xi[i].size() == 3);
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_EQ(s.isImpulseActive(i), impulse_status.isImpulseActive(i));
  }
  const Eigen::VectorXd lmd = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd gmm = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd beta = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd f_stack = Eigen::VectorXd::Random(impulse_status.dimp());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(impulse_status.dimp());
  const Eigen::VectorXd xi_stack = Eigen::VectorXd::Random(impulse_status.dimp());
  s.lmd = lmd;
  s.gmm = gmm;
  s.q = q;
  s.v = v;
  s.dv = dv;
  s.beta = beta;
  s.f_stack() = f_stack;
  s.mu_stack() = mu_stack;
  s.xi_stack() = xi_stack;
  EXPECT_TRUE(s.lmd.isApprox(lmd));
  EXPECT_TRUE(s.gmm.isApprox(gmm));
  EXPECT_TRUE(s.q.isApprox(q));
  EXPECT_TRUE(s.v.isApprox(v));
  EXPECT_TRUE(s.dv.isApprox(dv));
  EXPECT_TRUE(s.beta.isApprox(beta));
  EXPECT_TRUE(s.f_stack().isApprox(f_stack));
  EXPECT_TRUE(s.mu_stack().isApprox(mu_stack));
  EXPECT_TRUE(s.xi_stack().isApprox(xi_stack));
  s.set_f_vector();
  int dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_mu_vector();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_xi_vector();
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.xi[i].isApprox(s.xi_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s.f_stack().setZero();
  s.mu_stack().setZero();
  s.xi_stack().setZero();
  s.set_f_stack();
  s.set_mu_stack();
  s.set_xi_stack();
  dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.xi[i].isApprox(s.xi_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s = ImpulseSplitSolution(robot);
  s.setRandom(robot, impulse_status);
  EXPECT_TRUE(s.lmd.size() == robot.dimv());
  EXPECT_TRUE(s.gmm.size() == robot.dimv());
  EXPECT_TRUE(s.q.size() == robot.dimq());
  EXPECT_TRUE(s.v.size() == robot.dimv());
  EXPECT_TRUE(s.dv.size() == robot.dimv());
  EXPECT_TRUE(s.beta.size() == robot.dimv());
  EXPECT_TRUE(s.f.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.f[i].size() == 3);
  }
  EXPECT_TRUE(s.mu.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.mu[i].size() == 3);
  }
  EXPECT_TRUE(s.xi.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s.xi[i].size() == 3);
  }
  EXPECT_TRUE(s.f_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s.mu_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s.xi_stack().size() == impulse_status.dimp());
  EXPECT_EQ(s.dimf(), impulse_status.dimp());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_EQ(s.isImpulseActive(i), impulse_status.isImpulseActive(i));
  }
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.dv.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (impulse_status.hasActiveImpulse()) {
    EXPECT_FALSE(s.f_stack().isZero());
    EXPECT_FALSE(s.mu_stack().isZero());
    EXPECT_FALSE(s.xi_stack().isZero());
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_FALSE(s.f[i].isZero());
    }
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_FALSE(s.mu[i].isZero());
    }
  }
  const ImpulseSplitSolution s_random = ImpulseSplitSolution::Random(robot, impulse_status);
  EXPECT_TRUE(s_random.lmd.size() == robot.dimv());
  EXPECT_TRUE(s_random.gmm.size() == robot.dimv());
  EXPECT_TRUE(s_random.q.size() == robot.dimq());
  EXPECT_TRUE(s_random.v.size() == robot.dimv());
  EXPECT_TRUE(s_random.dv.size() == robot.dimv());
  EXPECT_TRUE(s_random.beta.size() == robot.dimv());
  EXPECT_TRUE(s_random.f.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s_random.f[i].size() == 3);
  }
  EXPECT_TRUE(s_random.mu.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s_random.mu[i].size() == 3);
  }
  EXPECT_TRUE(s_random.xi.size() == robot.max_point_contacts());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_TRUE(s_random.xi[i].size() == 3);
  }
  EXPECT_TRUE(s_random.f_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s_random.mu_stack().size() == impulse_status.dimp());
  EXPECT_TRUE(s_random.xi_stack().size() == impulse_status.dimp());
  EXPECT_EQ(s_random.dimf(), impulse_status.dimp());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_EQ(s_random.isImpulseActive(i), impulse_status.isImpulseActive(i));
  }
  EXPECT_FALSE(s_random.lmd.isZero());
  EXPECT_FALSE(s_random.gmm.isZero());
  EXPECT_FALSE(s_random.q.isZero());
  EXPECT_FALSE(s_random.v.isZero());
  EXPECT_FALSE(s_random.dv.isZero());
  EXPECT_FALSE(s_random.beta.isZero());
  if (impulse_status.hasActiveImpulse()) {
    EXPECT_FALSE(s_random.f_stack().isZero());
    EXPECT_FALSE(s_random.mu_stack().isZero());
    EXPECT_FALSE(s_random.xi_stack().isZero());
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s_random.isImpulseActive(i)) {
      EXPECT_FALSE(s_random.f[i].isZero());
    }
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (s_random.isImpulseActive(i)) {
      EXPECT_FALSE(s_random.mu[i].isZero());
    }
  }
}


void ImpulseSplitSolutionTest::TestIsApprox(const Robot& robot, 
                                            const ImpulseStatus& impulse_status) {
  ImpulseSplitSolution s(robot);
  s.setRandom(robot, impulse_status);
  EXPECT_FALSE(s.lmd.isZero());
  EXPECT_FALSE(s.gmm.isZero());
  EXPECT_FALSE(s.q.isZero());
  EXPECT_FALSE(s.v.isZero());
  EXPECT_FALSE(s.dv.isZero());
  EXPECT_FALSE(s.beta.isZero());
  if (impulse_status.hasActiveImpulse()) {
    EXPECT_FALSE(s.f_stack().isZero());
    EXPECT_FALSE(s.mu_stack().isZero());
    EXPECT_FALSE(s.xi_stack().isZero());
  }
  ImpulseSplitSolution s_ref = s;
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
  s_ref.dv.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.dv = s.dv;
  EXPECT_TRUE(s.isApprox(s_ref));
  s_ref.beta.setRandom();
  EXPECT_FALSE(s.isApprox(s_ref));
  s_ref.beta = s.beta;
  EXPECT_TRUE(s.isApprox(s_ref));
  if (impulse_status.hasActiveImpulse()) {
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
    s_ref.xi_stack().setRandom();
    EXPECT_FALSE(s.isApprox(s_ref));
    s_ref.xi_stack() = s.xi_stack();
    s_ref.set_xi_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
  else {
    s_ref.f_stack().setRandom();
    s_ref.set_f_vector();
    s_ref.mu_stack().setRandom();
    s_ref.set_mu_vector();
    s_ref.xi_stack().setRandom();
    s_ref.set_xi_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
}


void ImpulseSplitSolutionTest::TestIntegrate(const Robot& robot, 
                                             const ImpulseStatus& impulse_status) {
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitSolution s_ref = s;
  const double step_size = 0.3;
  s.integrate(robot, step_size, d);
  s_ref.lmd.noalias() += step_size * d.dlmd();
  s_ref.gmm.noalias() += step_size * d.dgmm();
  robot.integrateConfiguration(d.dq(), step_size, s_ref.q);
  s_ref.v.noalias() += step_size * d.dv();
  s_ref.dv.noalias() += step_size * d.ddv();
  s_ref.beta.noalias() += step_size * d.dbeta();
  if (impulse_status.hasActiveImpulse()) {
    s_ref.f_stack().noalias() += step_size * d.df();
    s_ref.set_f_vector();
    s_ref.mu_stack().noalias() += step_size * d.dmu();
    s_ref.set_mu_vector();
    s_ref.xi_stack().noalias() += step_size * d.dxi();
    s_ref.set_xi_vector();
  }
  EXPECT_TRUE(s.isApprox(s_ref));
}


TEST_F(ImpulseSplitSolutionTest, fixedBase) {
  Robot robot_without_impulse(fixed_base_urdf);
  TestWithoutImpulses(robot_without_impulse);
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::vector<bool> is_impulse_active = {false};
  ImpulseStatus impulse_status(is_impulse_active.size());
  impulse_status.setImpulseStatus(is_impulse_active);
  TestIsApprox(robot, impulse_status);
  TestIntegrate(robot, impulse_status);
  impulse_status.activateImpulse(0);
  TestWithImpulses(robot, impulse_status);
  TestIsApprox(robot, impulse_status);
  TestIntegrate(robot, impulse_status);
}


TEST_F(ImpulseSplitSolutionTest, floatingBase) {
  Robot robot_without_impulse(floating_base_urdf);
  TestWithoutImpulses(robot_without_impulse);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_impulse_active = {false, false, false, false};
  ImpulseStatus impulse_status(is_impulse_active.size());
  impulse_status.setImpulseStatus(is_impulse_active);
  TestIsApprox(robot, impulse_status);
  TestIntegrate(robot, impulse_status);
  std::random_device rnd;
  is_impulse_active.clear();
  for (const auto frame : contact_frames) {
    is_impulse_active.push_back(rnd()%2==0);
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  TestWithImpulses(robot, impulse_status);
  TestIsApprox(robot, impulse_status);
  TestIntegrate(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}