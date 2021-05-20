#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_switching_constraint_jacobian.hpp"
#include "idocp/ocp/split_switching_constraint_residual.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"
#include "idocp/riccati/backward_riccati_recursion_factorizer.hpp"
#include "idocp/riccati/split_constrained_riccati_factorization.hpp"
#include "idocp/riccati/split_riccati_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"


namespace idocp {

class SplitRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot);

  void testBackwardRecursion(const Robot& robot) const;
  void testBackwardRecursionWithSwitchingConstraint(const Robot& robot) const;
  void testForwardRecursion(const Robot& robot) const;

  double dt;
};


SplitRiccatiFactorization SplitRiccatiFactorizerTest::createRiccatiFactorization(const Robot& robot) {
  SplitRiccatiFactorization riccati(robot);
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  return riccati; 
}


void SplitRiccatiFactorizerTest::testBackwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitKKTMatrix kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  SplitKKTResidual kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.inverse();
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  Eigen::MatrixXd K(robot.dimu(), 2*robot.dimv());
  factorizer.getStateFeedbackGain(K);
  EXPECT_TRUE(lqr_policy_ref.K.isApprox(K));
  Eigen::MatrixXd Kq(dimu, dimv);
  Eigen::MatrixXd Kv(dimu, dimv);
  factorizer.getStateFeedbackGain(Kq, Kv);
  EXPECT_TRUE(lqr_policy_ref.K.leftCols(dimv).isApprox(Kq));
  EXPECT_TRUE(lqr_policy_ref.K.rightCols(dimv).isApprox(Kv));
}


void SplitRiccatiFactorizerTest::testBackwardRecursionWithSwitchingConstraint(const Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimi = impulse_status.dimf();
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitKKTMatrix kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  SplitKKTResidual kkt_residual = testhelper::CreateSplitKKTResidual(robot, impulse_status);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitStateConstraintJacobian jac(robot);
  jac.setImpulseStatus(impulse_status);
  jac.Phix().setRandom();
  jac.Phia().setRandom();
  jac.Phiu().setRandom();
  auto jac_ref = jac;
  SplitRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  SplitConstrainedRiccatiFactorization c_riccati(robot);
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, jac, riccati, c_riccati);

  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd GDtD = Eigen::MatrixXd::Zero(dimu+dimi, dimu+dimi);
  GDtD.topLeftCorner(dimu, dimu) = kkt_matrix.Quu;
  GDtD.topRightCorner(dimu, dimi) = jac_ref.Phiu().transpose();
  GDtD.bottomLeftCorner(dimi, dimu) = jac_ref.Phiu();
  const Eigen::MatrixXd GDtDinv = GDtD.inverse();
  Eigen::MatrixXd HtC = Eigen::MatrixXd::Zero(dimu+dimi, 2*dimv);
  HtC.topRows(dimu) = kkt_matrix_ref.Qxu.transpose();
  HtC.bottomRows(dimi) = jac_ref.Phix();
  Eigen::VectorXd htc = Eigen::VectorXd::Zero(dimu+dimi);
  htc.head(dimu) = kkt_residual_ref.lu;
  htc.tail(dimi) = kkt_residual_ref.P();
  const Eigen::MatrixXd KM = - GDtDinv * HtC;
  const Eigen::VectorXd km = - GDtDinv * htc;
  lqr_policy_ref.K = KM.topRows(dimu);
  lqr_policy_ref.k = km.head(dimu);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  Eigen::MatrixXd ODtD = GDtD;
  ODtD.topLeftCorner(dimu, dimu).setZero();
  const Eigen::MatrixXd MtKtODtDKM = KM.transpose() * ODtD * KM;
  riccati_ref.Pqq -= MtKtODtDKM.topLeftCorner(dimv, dimv);
  riccati_ref.Pqv -= MtKtODtDKM.topRightCorner(dimv, dimv);
  riccati_ref.Pvq -= MtKtODtDKM.bottomLeftCorner(dimv, dimv);
  riccati_ref.Pvv -= MtKtODtDKM.bottomRightCorner(dimv, dimv);
  riccati_ref.sq.noalias() -= jac_ref.Phix().transpose().topRows(dimv) * km.tail(dimi);
  riccati_ref.sv.noalias() -= jac_ref.Phix().transpose().bottomRows(dimv) * km.tail(dimi);
  constexpr double approx_eps = 1.0e-12;
  EXPECT_TRUE(c_riccati.DtM.isApprox((jac_ref.Phiu().transpose()*KM.bottomRows(dimi)), approx_eps));
  EXPECT_TRUE(c_riccati.KtDtM.isApprox((KM.topRows(dimu).transpose()*jac_ref.Phiu().transpose()*KM.bottomRows(dimi)), approx_eps));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq, approx_eps));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv, approx_eps));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq, approx_eps));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv, approx_eps));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq, approx_eps));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv, approx_eps));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  Eigen::MatrixXd K(robot.dimu(), 2*robot.dimv());
  factorizer.getStateFeedbackGain(K);
  EXPECT_TRUE(lqr_policy_ref.K.isApprox(K));
  Eigen::MatrixXd Kq(dimu, dimv);
  Eigen::MatrixXd Kv(dimu, dimv);
  factorizer.getStateFeedbackGain(Kq, Kv);
  EXPECT_TRUE(lqr_policy_ref.K.leftCols(dimv).isApprox(Kq, approx_eps));
  EXPECT_TRUE(lqr_policy_ref.K.rightCols(dimv).isApprox(Kv, approx_eps));
  EXPECT_TRUE(c_riccati.M().isApprox(KM.bottomRows(dimi), approx_eps));
  EXPECT_TRUE(c_riccati.m().isApprox(km.tail(dimi), approx_eps));
}


void SplitRiccatiFactorizerTest::testForwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_next_ref = riccati_next;
  SplitKKTMatrix kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  SplitKKTResidual kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  const Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.inverse();
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  SplitDirection d(robot);
  d.setRandom();
  SplitDirection d_next(robot);
  d_next.setRandom();
  SplitDirection d_ref = d;
  SplitDirection d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, dt, d, d_next);
  if (!robot.hasFloatingBase()) {
    kkt_matrix_ref.Fqq().setIdentity();
    kkt_matrix_ref.Fqv() = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  }
  d_ref.du = lqr_policy_ref.K * d_ref.dx + lqr_policy_ref.k;
  d_next_ref.dx = kkt_matrix_ref.Fxx * d.dx + kkt_residual_ref.Fx;
  d_next_ref.dv() += kkt_matrix_ref.Fvu * d_ref.du;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
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


TEST_F(SplitRiccatiFactorizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testBackwardRecursion(robot);
  testBackwardRecursionWithSwitchingConstraint(robot);
  testForwardRecursion(robot);
}


TEST_F(SplitRiccatiFactorizerTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testBackwardRecursion(robot);
  testBackwardRecursionWithSwitchingConstraint(robot);
  testForwardRecursion(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}