#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"
#include "idocp/riccati/unconstr_backward_riccati_recursion_factorizer.hpp"
#include "idocp/riccati/unconstr_riccati_factorizer.hpp"

#include "robot_factory.hpp"


namespace idocp {

class UnconstrRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();

    kkt_matrix = SplitKKTMatrix(robot);
    const Eigen::MatrixXd Qxx_seed = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
    kkt_matrix.Qxx = Qxx_seed * Qxx_seed.transpose();
    const Eigen::MatrixXd Quu_seed = Eigen::MatrixXd::Random(dimv, dimv);
    kkt_matrix.Qxu.setRandom();
    kkt_matrix.Quu = Quu_seed * Quu_seed.transpose();

    kkt_residual = SplitKKTResidual(robot);
    kkt_residual.lx.setRandom();
    kkt_residual.la.setRandom();
  }

  virtual void TearDown() {
  }

  SplitRiccatiFactorization createRiccatiFactorization() const;

  Robot robot;
  double dt;
  int dimv;
  SplitKKTMatrix kkt_matrix;
  SplitKKTResidual kkt_residual;
};


SplitRiccatiFactorization UnconstrRiccatiFactorizerTest::createRiccatiFactorization() const {
  SplitRiccatiFactorization riccati(robot);
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


TEST_F(UnconstrRiccatiFactorizerTest, backwardRiccatiRecursion) {
  auto riccati_next = createRiccatiFactorization();
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  UnconstrRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  UnconstrBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = createRiccatiFactorization();
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Qaa.inverse();
  lqr_policy_ref.K = - Ginv * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv * kkt_residual.la;
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
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(lqr_policy.isApprox(lqr_policy_ref));
}


TEST_F(UnconstrRiccatiFactorizerTest, forwardRiccatiRecursion) {
  auto riccati_next = createRiccatiFactorization();
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  UnconstrRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot);
  UnconstrBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = createRiccatiFactorization();
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  const Eigen::MatrixXd Ginv = kkt_matrix_ref.Quu.inverse();
  lqr_policy.K = - Ginv * kkt_matrix_ref.Qxu.transpose();
  lqr_policy.k = - Ginv * kkt_residual.la;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy, dt, riccati_ref);
  auto d = SplitDirection::Random(robot);
  auto d_next = SplitDirection::Random(robot);
  auto d_ref = d;
  auto d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(kkt_residual, dt, lqr_policy, d, d_next);
  d_ref.da() = lqr_policy.K * d_ref.dx + lqr_policy.k;
  Eigen::MatrixXd A(Eigen::MatrixXd::Zero(2*dimv, 2*dimv));
  A.topLeftCorner(dimv, dimv).setIdentity();
  A.topRightCorner(dimv, dimv).diagonal().fill(dt);
  A.bottomRightCorner(dimv, dimv).setIdentity();
  Eigen::MatrixXd B(Eigen::MatrixXd::Zero(2*dimv, dimv));
  B.bottomRows(dimv).diagonal().fill(dt);
  d_next_ref.dx = A * d_ref.dx + B * d_ref.da() + kkt_residual_ref.Fx;
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  factorizer.computeCostateDirection(riccati, d);
  Eigen::MatrixXd P(Eigen::MatrixXd::Zero(2*dimv, 2*dimv));
  P.topLeftCorner(dimv, dimv) = riccati.Pqq;
  P.topRightCorner(dimv, dimv) = riccati.Pqv;
  P.bottomLeftCorner(dimv, dimv) = riccati.Pqv.transpose();
  P.bottomRightCorner(dimv, dimv) = riccati.Pvv;
  d_ref.dlmdgmm = P * d_ref.dx;
  d_ref.dlmd() -= riccati.sq;
  d_ref.dgmm() -= riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}