#include <string>
#include <memory>
#include <limits>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"
#include "idocp/ocp/split_constrained_riccati_factorization.hpp"
#include "idocp/ocp/split_riccati_factorizer.hpp"

namespace idocp {

class SplitRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]) + std::sqrt(std::numeric_limits<double>::epsilon());
  }

  virtual void TearDown() {
  }

  SplitKKTMatrix createKKTMatrix(const Robot& robot) const;
  SplitKKTResidual createKKTResidual(const Robot& robot) const;
  SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot) const;

  void testBackwardRecursion(const Robot& robot) const;
  void testForwardRecursion(const Robot& robot) const;

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


SplitKKTMatrix SplitRiccatiFactorizerTest::createKKTMatrix(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  SplitKKTMatrix kkt_matrix(robot);
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
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
    kkt_matrix.Fqv().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  return kkt_matrix;
}


SplitKKTResidual SplitRiccatiFactorizerTest::createKKTResidual(const Robot& robot) const {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


SplitRiccatiFactorization SplitRiccatiFactorizerTest::createRiccatiFactorization(const Robot& robot) const {
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
  SplitKKTMatrix kkt_matrix = createKKTMatrix(robot);
  SplitKKTResidual kkt_residual = createKKTResidual(robot);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitRiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
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


void SplitRiccatiFactorizerTest::testForwardRecursion(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_next_ref = riccati_next;
  SplitKKTMatrix kkt_matrix = createKKTMatrix(robot);
  SplitKKTResidual kkt_residual = createKKTResidual(robot);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitRiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dtau, kkt_matrix_ref, kkt_residual_ref);
  const Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu().inverse();
  lqr_policy_ref.K = - Ginv  * kkt_matrix_ref.Qxu().transpose();
  lqr_policy_ref.k = - Ginv  * kkt_residual.lu();
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dtau, riccati_ref);
  SplitDirection d(robot);
  d.setRandom();
  SplitDirection d_next(robot);
  d_next.setRandom();
  SplitDirection d_ref = d;
  SplitDirection d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, dtau, d, d_next);
  if (!robot.hasFloatingBase()) {
    kkt_matrix_ref.Fqq().setIdentity();
  }
  d_ref.du() = lqr_policy_ref.K * d_ref.dx() + lqr_policy_ref.k;
  d_next_ref.dx() = kkt_matrix_ref.Fxx() * d.dx() + kkt_residual_ref.Fx();
  d_next_ref.dv() += kkt_matrix_ref.Fvu() * d_ref.du();
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(SplitRiccatiFactorizerTest, fixedBase) {
  Robot robot(fixed_base_urdf, {18});
  testBackwardRecursion(robot);
  testForwardRecursion(robot);
}


TEST_F(SplitRiccatiFactorizerTest, floating_base) {
  Robot robot(floating_base_urdf, {14, 24, 34, 44});
  testBackwardRecursion(robot);
  testForwardRecursion(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}