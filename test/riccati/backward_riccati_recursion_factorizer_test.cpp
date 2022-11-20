#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/backward_riccati_recursion_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace robotoc {

class BackwardRiccatiRecursionFactorizerTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }
  virtual void TearDown() {
  }
  double dt;
};


TEST_P(BackwardRiccatiRecursionFactorizerTest, test) {
  const auto robot = GetParam();
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
  bool has_next_sto_phase = true;
  factorizer.factorizeHamiltonian(riccati_next, kkt_matrix, riccati, has_next_sto_phase);
  riccati_ref.psi_x = kkt_matrix.hx + A.transpose() * riccati_next.P * kkt_matrix.fx + A.transpose() * riccati_next.Psi;
  riccati_ref.psi_u = kkt_matrix.hu + B.transpose() * riccati_next.P * kkt_matrix.fx + B.transpose() * riccati_next.Psi;
  riccati_ref.phi_x = A.transpose() * riccati_next.Phi;
  riccati_ref.phi_u = B.transpose() * riccati_next.Phi;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  has_next_sto_phase = false;
  factorizer.factorizeHamiltonian(riccati_next, kkt_matrix, riccati, has_next_sto_phase);
  riccati_ref.phi_x.setZero();
  riccati_ref.phi_u.setZero();
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  LQRPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, riccati);
  riccati_ref.P = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  riccati_ref.s = A.transpose() * riccati_next.s - A.transpose() * riccati_next.P * kkt_residual_ref.Fx - kkt_residual_ref.lx - H_ref * lqr_policy.k;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  has_next_sto_phase = true;
  factorizer.factorizeSTOFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, riccati, has_next_sto_phase);
  riccati_ref.Psi = riccati_ref.psi_x + lqr_policy.K.transpose() * riccati_ref.psi_u;
  riccati_ref.Phi = riccati_ref.phi_x + lqr_policy.K.transpose() * riccati_ref.phi_u;
  riccati_ref.xi = kkt_matrix.fx.transpose() * riccati_next.P * kkt_matrix.fx;
  riccati_ref.xi += kkt_matrix.Qtt;
  riccati_ref.xi += 2 * riccati_next.Psi.dot(kkt_matrix.fx);
  riccati_ref.xi += lqr_policy.T.dot(riccati_ref.psi_u);
  riccati_ref.xi += riccati_next.xi;
  riccati_ref.chi  = riccati_next.Phi.dot(kkt_matrix.fx);
  riccati_ref.chi += lqr_policy.T.dot(riccati_ref.phi_u);
  riccati_ref.chi += riccati_next.chi;
  riccati_ref.rho  = lqr_policy.W.dot(riccati_ref.phi_u);
  riccati_ref.rho += riccati_next.rho;
  riccati_ref.eta   = kkt_matrix.fx.transpose() * (riccati_next.P * kkt_residual.Fx - riccati_next.s);
  riccati_ref.eta  += kkt_residual.h;
  riccati_ref.eta  += riccati_next.Psi.dot(kkt_residual.Fx);
  riccati_ref.eta  += riccati_ref.psi_u.dot(lqr_policy.k);
  riccati_ref.eta  += riccati_next.eta;
  riccati_ref.iota  = riccati_next.Phi.dot(kkt_residual.Fx);
  riccati_ref.iota += riccati_ref.phi_u.dot(lqr_policy.k);
  riccati_ref.iota += riccati_next.iota;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  has_next_sto_phase = false;
  factorizer.factorizeSTOFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, riccati, has_next_sto_phase);
  riccati_ref.Phi.setZero();
  riccati_ref.chi  = 0.0;
  riccati_ref.rho  = 0.0;
  riccati_ref.iota = 0.0;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
}


TEST_P(BackwardRiccatiRecursionFactorizerTest, test_impact) {
  const auto robot = GetParam();
  const int dimv = robot.dimv();
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot);
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
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
  factorizer.factorizeSTOFactorization(riccati_next, kkt_matrix, kkt_residual, riccati);
  riccati_ref.Psi.setZero();
  riccati_ref.Phi  = A.transpose() * riccati_next.Phi;
  riccati_ref.xi   = 0.0;
  riccati_ref.chi  = 0.0;
  riccati_ref.rho  = riccati_next.rho;
  riccati_ref.eta  = 0.0;
  riccati_ref.iota = riccati_next.iota + riccati_next.Phi.dot(kkt_residual.Fx);
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, BackwardRiccatiRecursionFactorizerTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateQuadrupedalRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}