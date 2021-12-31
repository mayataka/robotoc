#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/unconstr_backward_riccati_recursion_factorizer.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace robotoc {

class UnconstrBackwardRiccatiRecursionFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateRobotManipulator();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();
    dimx = 2*robot.dimv();

    kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
    const Eigen::MatrixXd H_seed = Eigen::MatrixXd::Random(dimx+dimv, dimx+dimv);
    const Eigen::MatrixXd H = H_seed * H_seed.transpose();
    kkt_matrix.Qxx = H.topLeftCorner(dimx, dimx);
    kkt_matrix.Qxu = H.topRightCorner(dimx, dimv);
    kkt_matrix.Qaa = H.bottomRightCorner(dimv, dimv);
    kkt_matrix.Fxx.setZero();
    kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
  int dimv, dimx;
  SplitKKTMatrix kkt_matrix;
  SplitKKTResidual kkt_residual;
};


TEST_F(UnconstrBackwardRiccatiRecursionFactorizerTest, test) {
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
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
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx + A.transpose() * riccati_next.P * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu + A.transpose() * riccati_next.P * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Qaa + B.transpose() * riccati_next.P * B;
  const Eigen::VectorXd la_ref = B.transpose() * riccati_next.P * kkt_residual_ref.Fx - B.transpose() * riccati_next.s + kkt_residual_ref.la;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Qaa));
  EXPECT_TRUE(kkt_matrix.Qaa.isApprox(kkt_matrix.Qaa.transpose()));
  EXPECT_TRUE(la_ref.isApprox(kkt_residual.la));
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  LQRPolicy lqr_policy(robot);
  lqr_policy.K.setRandom();
  lqr_policy.k.setRandom();
  factorizer.factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, lqr_policy, dt, riccati);
  riccati_ref.P = F_ref - lqr_policy.K.transpose() * G_ref * lqr_policy.K;
  riccati_ref.s = A.transpose() * riccati_next.s - A.transpose() * riccati_next.P * kkt_residual_ref.Fx - kkt_residual_ref.lx - H_ref * lqr_policy.k;
  EXPECT_TRUE(riccati.isApprox(riccati_ref));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}