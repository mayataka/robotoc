#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_switching_constraint_jacobian.hpp"
#include "robotoc/ocp/split_switching_constraint_residual.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/backward_riccati_recursion_factorizer.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace robotoc {

class RiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test_backwardRecursion(const Robot& robot) const;
  void test_backwardRecursionWithSwitchingConstraint(const Robot& robot) const;
  void test_backwardRecursionImpulse(const Robot& robot) const;
  void test_forwardRecursion(const Robot& robot) const;
  void test_forwardRecursionImpulse(const Robot& robot) const;

  double dt;
};


void RiccatiFactorizerTest::test_backwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.llt().solve(Eigen::MatrixXd::Identity(dimu, dimu));
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(lqr_policy_ref.isApprox(lqr_policy));
}


void RiccatiFactorizerTest::test_backwardRecursionWithSwitchingConstraint(const Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimi = impulse_status.dimf();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  SplitSwitchingConstraintJacobian sc_jacobian(robot);
  SplitSwitchingConstraintResidual sc_residual(robot);
  sc_jacobian.setImpulseStatus(impulse_status);
  sc_residual.setImpulseStatus(impulse_status);
  sc_jacobian.Phix().setRandom();
  sc_jacobian.Phia().setRandom();
  sc_jacobian.Phiu().setRandom();
  sc_residual.P().setRandom();
  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  SplitConstrainedRiccatiFactorization c_riccati(robot);
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, 
                                      sc_jacobian, sc_residual, 
                                      riccati, c_riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd GDtD = Eigen::MatrixXd::Zero(dimu+dimi, dimu+dimi);
  GDtD.topLeftCorner(dimu, dimu) = kkt_matrix.Quu;
  GDtD.topRightCorner(dimu, dimi) = sc_jacobian.Phiu().transpose();
  GDtD.bottomLeftCorner(dimi, dimu) = sc_jacobian.Phiu();
  const Eigen::MatrixXd GDtDinv = GDtD.inverse();
  Eigen::MatrixXd HtC = Eigen::MatrixXd::Zero(dimu+dimi, 2*dimv);
  HtC.topRows(dimu) = kkt_matrix_ref.Qxu.transpose();
  HtC.bottomRows(dimi) = sc_jacobian.Phix();
  Eigen::VectorXd htc = Eigen::VectorXd::Zero(dimu+dimi);
  htc.head(dimu) = kkt_residual_ref.lu;
  htc.tail(dimi) = sc_residual.P();
  const Eigen::MatrixXd KM = - GDtDinv * HtC;
  const Eigen::VectorXd km = - GDtDinv * htc;
  lqr_policy_ref.K = KM.topRows(dimu);
  lqr_policy_ref.k = km.head(dimu);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  Eigen::MatrixXd ODtD = GDtD;
  ODtD.topLeftCorner(dimu, dimu).setZero();
  const Eigen::MatrixXd MtKtODtDKM = KM.transpose() * ODtD * KM;
  riccati_ref.P -= MtKtODtDKM;
  riccati_ref.s -= sc_jacobian.Phix().transpose() * km.tail(dimi);
  EXPECT_TRUE(c_riccati.DtM.isApprox((sc_jacobian.Phiu().transpose()*KM.bottomRows(dimi))));
  EXPECT_TRUE(c_riccati.KtDtM.isApprox((KM.topRows(dimu).transpose()*sc_jacobian.Phiu().transpose()*KM.bottomRows(dimi))));
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(c_riccati.M().isApprox(KM.bottomRows(dimi)));
  EXPECT_TRUE(c_riccati.m().isApprox(km.tail(dimi)));
  EXPECT_TRUE(lqr_policy_ref.isApprox(lqr_policy));
}


void RiccatiFactorizerTest::test_backwardRecursionImpulse(const Robot& robot) const {
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateImpulseSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateImpulseSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
}


void RiccatiFactorizerTest::test_forwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  const Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.llt().solve(Eigen::MatrixXd::Identity(dimu, dimu));
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, riccati_ref);
  auto d = SplitDirection::Random(robot);
  auto d_next = SplitDirection::Random(robot);
  auto d_ref = d;
  auto d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, lqr_policy, d, d_next);
  d_ref.du = lqr_policy_ref.K * d_ref.dx + lqr_policy_ref.k;
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx;
  d_next_ref.dv() += kkt_matrix_ref.Fvu * d_ref.du;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmdgmm = riccati.P * d.dx - riccati.s;
  EXPECT_TRUE(d.isApprox(d_ref));
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  SplitConstrainedRiccatiFactorization c_riccati(robot);
  c_riccati.setImpulseStatus(impulse_status.dimf());
  c_riccati.M().setRandom();
  c_riccati.m().setRandom();
  d.setImpulseStatus(impulse_status);
  d.dxi().setRandom();
  d_ref = d;
  factorizer.computeLagrangeMultiplierDirection(c_riccati, d);
  d_ref.dxi() = c_riccati.M() * d.dx + c_riccati.m();
  EXPECT_TRUE(d.isApprox(d_ref));
}


void RiccatiFactorizerTest::test_forwardRecursionImpulse(const Robot& robot) const {
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix = testhelper::CreateImpulseSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateImpulseSplitKKTResidual(robot);
  RiccatiFactorizer factorizer(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  auto d = ImpulseSplitDirection::Random(robot);
  auto d_next = SplitDirection::Random(robot);
  auto d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d, d_next);
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx;
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  auto d_ref = d;
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmdgmm = riccati.P * d.dx - riccati.s;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(RiccatiFactorizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test_backwardRecursion(robot);
  test_backwardRecursionWithSwitchingConstraint(robot);
  test_backwardRecursionImpulse(robot);
  test_forwardRecursion(robot);
  test_forwardRecursionImpulse(robot);
}


TEST_F(RiccatiFactorizerTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test_backwardRecursion(robot);
  test_backwardRecursionWithSwitchingConstraint(robot);
  test_backwardRecursionImpulse(robot);
  test_forwardRecursion(robot);
  test_forwardRecursionImpulse(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}