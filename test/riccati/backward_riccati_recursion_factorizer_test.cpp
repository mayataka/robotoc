#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"
#include "idocp/riccati/backward_riccati_recursion_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace idocp {

class BackwardRiccatiRecursionFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  void test_impulse(const Robot& robot) const;

  double dt;
};


void BackwardRiccatiRecursionFactorizerTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  const auto kkt_matrix_ref = kkt_matrix;
  const auto kkt_residual_ref = kkt_residual;
  BackwardRiccatiRecursionFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  const Eigen::MatrixXd A = kkt_matrix_ref.Fxx;
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimu);
  B.bottomRows(dimv) = kkt_matrix_ref.Fvu;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * riccati_next.P * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu + A.transpose() * riccati_next.P * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Quu + B.transpose() * riccati_next.P * B;
  const Eigen::VectorXd lu_ref = B.transpose() * riccati_next.P * kkt_residual_ref.Fx - B.transpose() * riccati_next.s + kkt_residual_ref.lu;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Quu));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));
  EXPECT_TRUE(lu_ref.isApprox(kkt_residual.lu));
  SplitRiccatiFactorization riccati(robot), riccati_ref(robot);
  LQRPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeHamiltonian(kkt_matrix, riccati);
  const Eigen::VectorXd hx_ref = kkt_matrix.hx + A.transpose() * riccati_next.P * kkt_matrix.fx;
  const Eigen::VectorXd hu_ref = kkt_matrix.hu + B.transpose() * riccati_next.P * kkt_matrix.fx;
  EXPECT_TRUE(kkt_matrix.hx.isApprox(hx_ref));
  EXPECT_TRUE(kkt_matrix.hu.isApprox(hu_ref));
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, riccati);
  riccati_ref.P = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  riccati_ref.s = A.transpose() * riccati_next.s - A.transpose() * riccati_next.P * kkt_residual_ref.Fx - kkt_residual_ref.lx - H_ref * lqr_policy.k;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
}


void BackwardRiccatiRecursionFactorizerTest::test_impulse(const Robot& robot) const {
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateImpulseSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateImpulseSplitKKTResidual(robot);
  const auto kkt_matrix_ref = kkt_matrix;
  const auto kkt_residual_ref = kkt_residual;
  BackwardRiccatiRecursionFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, kkt_matrix);
  const Eigen::MatrixXd A = kkt_matrix.Fxx;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * riccati_next.P * A;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  SplitRiccatiFactorization riccati(robot), riccati_ref(robot);
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, riccati);
  riccati_ref.P = F_ref;
  riccati_ref.s = A.transpose() * riccati_next.s - A.transpose() * riccati_next.P * kkt_residual_ref.Fx - kkt_residual_ref.lx;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
}


TEST_F(BackwardRiccatiRecursionFactorizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
  test_impulse(robot);
}


TEST_F(BackwardRiccatiRecursionFactorizerTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
  test_impulse(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}