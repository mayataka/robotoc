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

#include "robot_factory.hpp"


namespace idocp {

class UnconstrBackwardRiccatiRecursionFactorizerTest : public ::testing::Test {
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


SplitRiccatiFactorization UnconstrBackwardRiccatiRecursionFactorizerTest::createRiccatiFactorization() const {
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


TEST_F(UnconstrBackwardRiccatiRecursionFactorizerTest, test) {
  const auto riccati_next = createRiccatiFactorization();
  const auto kkt_matrix_ref = kkt_matrix;
  const auto kkt_residual_ref = kkt_residual;
  UnconstrBackwardRiccatiRecursionFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, dt, kkt_matrix, kkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv).setIdentity();
  A.topRightCorner(dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  A.bottomLeftCorner(dimv, dimv).setZero();
  A.bottomRightCorner(dimv, dimv).setIdentity();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimv);
  B.topRows(dimv).setZero();
  B.bottomRows(dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * P_next * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu + A.transpose() * P_next * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Qaa + B.transpose() * P_next * B;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  const Eigen::VectorXd la_ref = B.transpose() * P_next * kkt_residual_ref.Fx - B.transpose() * sx_next + kkt_residual_ref.la;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Qaa));
  EXPECT_TRUE(kkt_matrix.Qaa.isApprox(kkt_matrix.Qaa.transpose()));
  EXPECT_TRUE(la_ref.isApprox(kkt_residual.la));
  auto riccati = createRiccatiFactorization();
  LQRPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, dt, riccati);
  const Eigen::MatrixXd P_ref = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * kkt_residual_ref.Fx - kkt_residual_ref.lx - H_ref * lqr_policy.k;
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}