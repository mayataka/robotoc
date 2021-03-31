#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/unocp/backward_unriccati_recursion_factorizer.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SplitUnRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();

    unkkt_matrix = SplitUnKKTMatrix(robot);
    const Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
    unkkt_matrix.Q = seed * seed.transpose();

    unkkt_residual = SplitUnKKTResidual(robot);
    unkkt_residual.KKT_residual.setRandom();
  }

  virtual void TearDown() {
  }

  SplitRiccatiFactorization createRiccatiFactorization() const;

  Robot robot;
  double dt;
  int dimv;
  SplitUnKKTMatrix unkkt_matrix;
  SplitUnKKTResidual unkkt_residual;
};


SplitRiccatiFactorization SplitUnRiccatiFactorizerTest::createRiccatiFactorization() const {
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


TEST_F(SplitUnRiccatiFactorizerTest, backwardRiccatiRecursion) {
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization();
  SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  SplitUnRiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardUnRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization();
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, unkkt_matrix, unkkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, unkkt_matrix_ref, unkkt_residual_ref);
  Eigen::MatrixXd Ginv = unkkt_matrix_ref.Qaa().inverse();
  lqr_policy_ref.K = - Ginv * unkkt_matrix_ref.Qax();
  lqr_policy_ref.k = - Ginv * unkkt_residual.la();
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, unkkt_matrix_ref, unkkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(unkkt_matrix.Qxx().isApprox(unkkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(unkkt_matrix.Qaa().isApprox(unkkt_matrix.Qaa().transpose()));
  Eigen::MatrixXd K(robot.dimv(), 2*robot.dimv());
  factorizer.getStateFeedbackGain(K);
  EXPECT_TRUE(lqr_policy_ref.K.isApprox(K));
  Eigen::MatrixXd Kq(dimv, dimv);
  Eigen::MatrixXd Kv(dimv, dimv);
  factorizer.getStateFeedbackGain(Kq, Kv);
  EXPECT_TRUE(lqr_policy_ref.K.leftCols(dimv).isApprox(Kq));
  EXPECT_TRUE(lqr_policy_ref.K.rightCols(dimv).isApprox(Kv));
}


TEST_F(SplitUnRiccatiFactorizerTest, forwardRiccatiRecursion) {
  SplitRiccatiFactorization riccati_next = createRiccatiFactorization();
  SplitRiccatiFactorization riccati_next_ref = riccati_next;
  SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  SplitUnRiccatiFactorizer factorizer(robot);
  LQRStateFeedbackPolicy lqr_policy_ref(robot);
  BackwardUnRiccatiRecursionFactorizer backward_recursion_ref(robot);
  SplitRiccatiFactorization riccati = createRiccatiFactorization();
  SplitRiccatiFactorization riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, unkkt_matrix, unkkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, unkkt_matrix_ref, unkkt_residual_ref);
  const Eigen::MatrixXd Ginv = unkkt_matrix_ref.Qaa().inverse();
  lqr_policy_ref.K = - Ginv * unkkt_matrix_ref.Qax();
  lqr_policy_ref.k = - Ginv * unkkt_residual.la();
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, unkkt_matrix_ref, unkkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  SplitDirection d(robot);
  d.setRandom();
  SplitDirection d_next(robot);
  d_next.setRandom();
  SplitDirection d_ref = d;
  SplitDirection d_next_ref = d_next;
  factorizer.forwardRiccatiRecursion(unkkt_residual, d, dt, d_next);
  d_ref.da() = lqr_policy_ref.K * d_ref.dx() + lqr_policy_ref.k;
  Eigen::MatrixXd A(Eigen::MatrixXd::Zero(2*dimv, 2*dimv));
  A.topLeftCorner(dimv, dimv).setIdentity();
  A.topRightCorner(dimv, dimv).diagonal().fill(dt);
  A.bottomRightCorner(dimv, dimv).setIdentity();
  Eigen::MatrixXd B(Eigen::MatrixXd::Zero(2*dimv, dimv));
  B.bottomRows(dimv).diagonal().fill(dt);
  d_next_ref.dx() = A * d_ref.dx() + B * d_ref.da() + unkkt_residual_ref.Fx();
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  factorizer.computeCostateDirection(riccati, d);
  Eigen::MatrixXd P(Eigen::MatrixXd::Zero(2*dimv, 2*dimv));
  P.topLeftCorner(dimv, dimv) = riccati.Pqq;
  P.topRightCorner(dimv, dimv) = riccati.Pqv;
  P.bottomLeftCorner(dimv, dimv) = riccati.Pqv.transpose();
  P.bottomRightCorner(dimv, dimv) = riccati.Pvv;
  d_ref.dlmdgmm() = P * d_ref.dx();
  d_ref.dlmd() -= riccati.sq;
  d_ref.dgmm() -= riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}