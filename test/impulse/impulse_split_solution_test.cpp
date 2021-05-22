#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/ocp/split_solution.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseSplitSolutionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIntegrate(const Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseSplitSolutionTest::test(const Robot& robot, const ImpulseStatus& impulse_status) { 
  std::random_device rnd;
  ImpulseSplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.dv.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxPointContacts());
  EXPECT_EQ(s.mu.size(), robot.maxPointContacts());
  EXPECT_EQ(s.f_stack().size(), impulse_status.dimf());
  EXPECT_EQ(s.mu_stack().size(), impulse_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.f[i].size(), 3);
    EXPECT_EQ(s.mu[i].size(), 3);
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.isImpulseActive(i), impulse_status.isImpulseActive(i));
  }
  EXPECT_TRUE(s.lmd.isZero());
  EXPECT_TRUE(s.gmm.isZero());
  if (robot.hasFloatingBase()) {
    Eigen::VectorXd q_ref = Eigen::VectorXd::Zero(robot.dimq());
    q_ref(6) = 1.0;
    EXPECT_TRUE(s.q.isApprox(q_ref));
    std::cout << s.q.transpose() << std::endl;
  }
  else {
    EXPECT_TRUE(s.q.isZero());
  }
  EXPECT_TRUE(s.v.isZero());
  EXPECT_TRUE(s.dv.isZero());
  EXPECT_TRUE(s.beta.isZero());
  EXPECT_TRUE(s.f_stack().isZero());
  EXPECT_TRUE(s.mu_stack().isZero());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_TRUE(s.f[i].isZero());
    EXPECT_TRUE(s.mu[i].isZero());
  }
  const Eigen::VectorXd f_stack = Eigen::VectorXd::Random(impulse_status.dimf());
  const Eigen::VectorXd mu_stack = Eigen::VectorXd::Random(impulse_status.dimf());
  s.f_stack() = f_stack;
  s.mu_stack() = mu_stack;
  s.set_f_vector();
  int dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  s.set_mu_vector();
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
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
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.f[i].isApprox(s.f_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    if (s.isImpulseActive(i)) {
      EXPECT_TRUE(s.mu[i].isApprox(s.mu_stack().segment<3>(dimf_stack)));
      dimf_stack += 3;
    }
  }
  s = ImpulseSplitSolution::Random(robot, impulse_status);
  EXPECT_EQ(s.lmd.size(), robot.dimv());
  EXPECT_EQ(s.gmm.size(), robot.dimv());
  EXPECT_EQ(s.q.size(), robot.dimq());
  EXPECT_EQ(s.v.size(), robot.dimv());
  EXPECT_EQ(s.dv.size(), robot.dimv());
  EXPECT_EQ(s.beta.size(), robot.dimv());
  EXPECT_EQ(s.f.size(), robot.maxPointContacts());
  EXPECT_EQ(s.mu.size(), robot.maxPointContacts());
  EXPECT_EQ(s.f_stack().size(), impulse_status.dimf());
  EXPECT_EQ(s.mu_stack().size(), impulse_status.dimf());
  EXPECT_EQ(s.dimi(), impulse_status.dimf());
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    EXPECT_EQ(s.f[i].size(), 3);
    EXPECT_EQ(s.mu[i].size(), 3);
  }
  for (int i=0; i<robot.maxPointContacts(); ++i) {
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
  }
}


void ImpulseSplitSolutionTest::testIsApprox(const Robot& robot, 
                                            const ImpulseStatus& impulse_status) {
  auto s = ImpulseSplitSolution::Random(robot, impulse_status);
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
  }
  else {
    s_ref.f_stack().setRandom();
    s_ref.set_f_vector();
    s_ref.mu_stack().setRandom();
    s_ref.set_mu_vector();
    EXPECT_TRUE(s.isApprox(s_ref));
  }
}


void ImpulseSplitSolutionTest::testIntegrate(const Robot& robot, 
                                             const ImpulseStatus& impulse_status) {
  auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  const auto d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto s_ref = s;
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
  }
  EXPECT_TRUE(s.isApprox(s_ref));
}


TEST_F(ImpulseSplitSolutionTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  testIntegrate(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  testIntegrate(robot, impulse_status);
}


TEST_F(ImpulseSplitSolutionTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  testIntegrate(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  testIntegrate(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}