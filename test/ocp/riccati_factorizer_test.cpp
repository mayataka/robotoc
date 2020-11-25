#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"

namespace idocp {

class RiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  KKTMatrix createKKTMatrix(const Robot& robot) const;
  KKTResidual createKKTResidual(const Robot& robot) const;
  RiccatiFactorization createRiccatiFactorization(const Robot& robot) const;

  void testBackwardRecursion(const Robot& robot) const;
  void testForwardRecursionWithoutStateConstraint(const Robot& robot) const;
  void testForwardRecursionWithStateConstraint(const Robot& robot) const;
  void testFactorizeStateConstraintFactorization(const Robot& robot) const;

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


KKTMatrix RiccatiFactorizerTest::createKKTMatrix(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq() = seed * seed.transpose();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Qvv() = seed * seed.transpose();
  kkt_matrix.Qqu().setRandom();
  kkt_matrix.Qvu().setRandom();
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  kkt_matrix.Quu() = seed * seed.transpose();
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.has_floating_base()) {
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  return kkt_matrix;
}


KKTResidual RiccatiFactorizerTest::createKKTResidual(const Robot& robot) const {
  KKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


RiccatiFactorization RiccatiFactorizerTest::createRiccatiFactorization(const Robot& robot) const {
  RiccatiFactorization riccati(robot);
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  riccati.Pi.setRandom();
  riccati.pi.setRandom();
  riccati.N.setRandom();
  riccati.N = (riccati.N + riccati.N.transpose()).eval();
  riccati.n.setRandom();
  return riccati; 
}


void RiccatiFactorizerTest::testBackwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  KKTMatrix kkt_matrix = createKKTMatrix(robot);
  KKTResidual kkt_residual = createKKTResidual(robot);
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  KKTResidual kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  RiccatiFactorization riccati = createRiccatiFactorization(robot);
  RiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu().inverse();
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu().transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu();
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dtau, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
  Eigen::MatrixXd K(robot.dimu(), 2*robot.dimv());
  factorizer.getStateFeedbackGain(K);
  EXPECT_TRUE(lqr_policy_ref.K.isApprox(K));
  Eigen::MatrixXd Kq(dimu, dimv);
  Eigen::MatrixXd Kv(dimu, dimv);
  factorizer.getStateFeedbackGain(Kq, Kv);
  EXPECT_TRUE(lqr_policy_ref.K.leftCols(dimv).isApprox(Kq));
  EXPECT_TRUE(lqr_policy_ref.K.rightCols(dimv).isApprox(Kv));
}


void RiccatiFactorizerTest::testForwardRecursionWithoutStateConstraint(const Robot& robot) const {
  // const int dimv = robot.dimv();
  // const int dimu = robot.dimu();
  // const RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  // KKTMatrix kkt_matrix = createKKTMatrix(robot);
  // KKTResidual kkt_residual = createKKTResidual(robot);
  // KKTMatrix kkt_matrix_ref = kkt_matrix;
  // KKTResidual kkt_residual_ref = kkt_residual;
  // const SplitDirection d = SplitDirection::Random(robot);
  // SplitDirection d_next(robot), d_next_ref(robot);
  // RiccatiFactorizer factorizer(robot);
  // factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d, dtau, d_next);
  // d_next_ref.dx() = kkt_matrix.Fxx() * d.dx() + kkt_matrix.Fxu() * d.du() + kkt_residual.Fx();
  // EXPECT_TRUE(d_next.isApprox(d_next_ref));
}


void RiccatiFactorizerTest::testForwardRecursionWithStateConstraint(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  RiccatiFactorization riccati_next_ref = riccati_next;
  KKTMatrix kkt_matrix = createKKTMatrix(robot);
  KKTResidual kkt_residual = createKKTResidual(robot);
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  KKTResidual kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  RiccatiFactorization riccati = createRiccatiFactorization(robot);
  RiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu().inverse();
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu().transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu();
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dtau, riccati_ref);
  factorizer.forwardRiccatiRecursionParallel(kkt_matrix, kkt_residual, true);
  kkt_matrix_ref.Fxx() += kkt_matrix_ref.Fxu() * lqr_policy_ref.K;
  kkt_residual_ref.Fx() += kkt_matrix_ref.Fxu() * lqr_policy_ref.k;
  const Eigen::MatrixXd BGinvBt = kkt_matrix_ref.Fxu() * Ginv * kkt_matrix_ref.Fxu().transpose();
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  factorizer.forwardRiccatiRecursionSerial(riccati, kkt_matrix, kkt_residual, dtau, riccati_next, true);
  riccati_next_ref.Pi = kkt_matrix_ref.Fxx() * riccati.Pi;
  riccati_next_ref.pi = kkt_residual_ref.Fx() + kkt_matrix_ref.Fxx() * riccati.pi;
  riccati_next_ref.N = BGinvBt + kkt_matrix_ref.Fxx() * riccati.N * kkt_matrix_ref.Fxx().transpose();
  EXPECT_TRUE(riccati_next.Pi.isApprox(riccati_next_ref.Pi));
  EXPECT_TRUE(riccati_next.pi.isApprox(riccati_next_ref.pi));
  EXPECT_TRUE(riccati_next.N.isApprox(riccati_next_ref.N));
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*dimv);
  SplitDirection d(robot);
  d.setRandom();
  SplitDirection d_ref = d;
  factorizer.computeStateDirection(riccati, dx0, d, false);
  d_ref.dx().noalias() = riccati.Pi * dx0 + riccati.pi;
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeCostateDirection(riccati, d, false);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeControlInputDirection(riccati, d, false);
  d_ref.du() = lqr_policy_ref.K * d.dx() + lqr_policy_ref.k;
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeStateDirection(riccati, dx0, d, true);
  d_ref.dx().noalias() = riccati.Pi * dx0 + riccati.pi - riccati.N * riccati.n;
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeCostateDirection(riccati, d, true);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq + riccati.n.head(dimv);
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv + riccati.n.tail(dimv);
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeControlInputDirection(riccati, d, true);
  d_ref.du() = lqr_policy_ref.K * d.dx() + lqr_policy_ref.k - Ginv * kkt_matrix_ref.Fxu().transpose() * riccati.n;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void RiccatiFactorizerTest::testFactorizeStateConstraintFactorization(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const KKTMatrix kkt_matrix = createKKTMatrix(robot);
  const KKTResidual kkt_residual = createKKTResidual(robot);
  RiccatiFactorizer factorizer(robot);
  const int dimf = 6;
  Eigen::MatrixXd T_next(2*dimv, dimf), T(2*dimv, dimf), T_ref(2*dimv, dimf);
  T_next.setRandom();
  T.setRandom();
  T_ref = T;
  factorizer.backwardStateConstraintFactorization(T_next, kkt_matrix, dtau, T);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  if (robot.has_floating_base()) {
    A.topLeftCorner(dimv, dimv) = kkt_matrix.Fqq();
  }
  else {
    A.topLeftCorner(dimv, dimv).setIdentity();
  }
  A.topRightCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv();
  T_ref = A.transpose() * T_next;
  EXPECT_TRUE(T.isApprox(T_ref));
}


TEST_F(RiccatiFactorizerTest, fixedBase) {
  testBackwardRecursion(fixed_base_robot);
  testForwardRecursionWithoutStateConstraint(fixed_base_robot);
  testForwardRecursionWithStateConstraint(fixed_base_robot);
  testFactorizeStateConstraintFactorization(fixed_base_robot);
}


TEST_F(RiccatiFactorizerTest, floating_base) {
  testBackwardRecursion(floating_base_robot);
  testForwardRecursionWithoutStateConstraint(floating_base_robot);
  testForwardRecursionWithStateConstraint(floating_base_robot);
  testFactorizeStateConstraintFactorization(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}