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
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace idocp {

class UnconstrRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
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


TEST_F(UnconstrRiccatiFactorizerTest, backwardRiccatiRecursion) {
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  UnconstrRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot), lqr_policy_ref(robot);
  UnconstrBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_ref = riccati;
  factorizer.backwardRiccatiRecursion(riccati_next, dt, kkt_matrix, kkt_residual, riccati, lqr_policy);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, dt, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd Ginv = kkt_matrix_ref.Qaa.inverse();
  lqr_policy_ref.K = - Ginv * kkt_matrix_ref.Qxu.transpose();
  lqr_policy_ref.k = - Ginv * kkt_residual.la;
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, lqr_policy_ref, dt, riccati_ref);
  EXPECT_TRUE(riccati.P.isApprox(riccati_ref.P));
  EXPECT_TRUE(riccati.s.isApprox(riccati_ref.s));
  EXPECT_TRUE(riccati.P.isApprox(riccati.P.transpose()));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(lqr_policy.isApprox(lqr_policy_ref));
}


TEST_F(UnconstrRiccatiFactorizerTest, forwardRiccatiRecursion) {
  const auto riccati_next = testhelper::CreateSplitRiccatiFactorization(robot);
  auto riccati_next_ref = riccati_next;
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  UnconstrRiccatiFactorizer factorizer(robot);
  LQRPolicy lqr_policy(robot);
  UnconstrBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  auto riccati = testhelper::CreateSplitRiccatiFactorization(robot);
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
  d_ref.dlmdgmm = riccati.P * d_ref.dx;
  d_ref.dlmdgmm -= riccati.s;
  EXPECT_TRUE(d.isApprox(d_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}