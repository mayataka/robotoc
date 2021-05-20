#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/riccait/split_riccati_factorization.hpp"
#include "idocp/riccait/lqr_policy.hpp"
#include "idocp/riccait/backward_riccati_recursion_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"


namespace idocp {

class BackwardRiccatiRecursionFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  SplitRiccatiFactorization createRiccatiFactorization(const Robot& robot) const;

  void test(const Robot& robot) const;

  double dt;
};


SplitRiccatiFactorization BackwardRiccatiRecursionFactorizerTest::createRiccatiFactorization(const Robot& robot) const {
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


void BackwardRiccatiRecursionFactorizerTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const SplitRiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  SplitKKTMatrix kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  SplitKKTResidual kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  const SplitKKTResidual kkt_residual_ref = kkt_residual;
  BackwardRiccatiRecursionFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, dt, kkt_matrix, kkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv).setIdentity();
  if (robot.hasFloatingBase()) {
    A.topLeftCorner(robot.dim_passive(), robot.dim_passive()) 
        = kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive());
  }
  A.topRightCorner(dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    A.topRightCorner(dimv, dimv).topLeftCorner(robot.dim_passive(), robot.dim_passive()) 
        = kkt_matrix.Fqv().topLeftCorner(robot.dim_passive(), robot.dim_passive());
  }
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimu);
  B.bottomRows(dimv) = kkt_matrix_ref.Fvu;
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * P_next * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu + A.transpose() * P_next * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Quu + B.transpose() * P_next * B;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  const Eigen::VectorXd lu_ref = B.transpose() * P_next * kkt_residual_ref.Fx - B.transpose() * sx_next + kkt_residual_ref.lu;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Quu));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(lu_ref.isApprox(kkt_residual.lu));
  SplitRiccatiFactorization riccati = createRiccatiFactorization(robot);
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


TEST_F(BackwardRiccatiRecursionFactorizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(BackwardRiccatiRecursionFactorizerTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}