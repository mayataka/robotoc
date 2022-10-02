#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/unconstr/unconstr_ocp.hpp"
#include "robotoc/riccati/unconstr_riccati_factorizer.hpp"
#include "robotoc/riccati/unconstr_riccati_recursion.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"
#include "direction_factory.hpp"
#include "kkt_factory.hpp"
#include "riccati_factory.hpp"


namespace robotoc {

class UnconstrRiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateRobotManipulator();
    dimv = robot.dimv();
    dimx = 2 * robot.dimv();
    N = 20;
    T = 1;
    dt = T / N;

    kkt_matrix = KKTMatrix(N+1, SplitKKTMatrix(robot));
    kkt_residual = KKTResidual(N+1, SplitKKTResidual(robot));
    d = Direction(N+1, SplitDirection(robot));
    riccati_factorization = UnconstrRiccatiFactorization(N+1, SplitRiccatiFactorization(robot));
    for (int i=0; i<=N; ++i) {
      const Eigen::MatrixXd H_seed = Eigen::MatrixXd::Random(dimx+dimv, dimx+dimv);
      const Eigen::MatrixXd H = H_seed * H_seed.transpose();
      kkt_matrix[i].Qxx = H.topLeftCorner(dimx, dimx);
      kkt_matrix[i].Qxu = H.topRightCorner(dimx, dimv);
      kkt_matrix[i].Qaa = H.bottomRightCorner(dimv, dimv);
      kkt_residual[i] = testhelper::CreateSplitKKTResidual(robot);
      d[i].setRandom();
      riccati_factorization[i] = testhelper::CreateSplitRiccatiFactorization(robot);
    }
  }

  virtual void TearDown() {
  }

  Robot robot;
  int N, dimv, dimx;
  double T, dt;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  Direction d;
  UnconstrRiccatiFactorization riccati_factorization;
};


TEST_F(UnconstrRiccatiRecursionTest, test) {
  auto riccati_factorization_ref = riccati_factorization;
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto ocp = OCP(robot, cost, constraints, T, N);
  UnconstrRiccatiRecursion riccati_recursion(ocp);
  riccati_recursion.backwardRiccatiRecursion(kkt_matrix, kkt_residual, 
                                             riccati_factorization);
  riccati_factorization_ref[N].P = kkt_matrix_ref[N].Qxx;
  riccati_factorization_ref[N].s = - kkt_residual_ref[N].lx;
  EXPECT_TRUE(riccati_factorization[N].isApprox(riccati_factorization_ref[N]));
  UnconstrRiccatiFactorizer factorizer(robot);
  std::vector<LQRPolicy> lqr_policy(N, LQRPolicy(robot));
  for (int i=N-1; i>=0; --i) {
    factorizer.backwardRiccatiRecursion(
        riccati_factorization_ref[i+1], dt, kkt_matrix_ref[i], kkt_residual_ref[i], 
        riccati_factorization_ref[i], lqr_policy[i]);
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(kkt_matrix[i].isApprox(kkt_matrix_ref[i]));
    EXPECT_TRUE(kkt_residual[i].isApprox(kkt_residual_ref[i]));
    EXPECT_TRUE(riccati_factorization[i].isApprox(riccati_factorization_ref[i]));
  }
  d[0].dx.setRandom();
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(kkt_residual, riccati_factorization, d);
  for (int i=0; i<N; ++i) {
    factorizer.forwardRiccatiRecursion(kkt_residual_ref[i], dt, lqr_policy[i],
                                       d_ref[i], d_ref[i+1]);
    factorizer.computeCostateDirection(riccati_factorization_ref[i], d_ref[i]);
  }
  factorizer.computeCostateDirection(riccati_factorization_ref[N], d_ref[N]);
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(d[i].isApprox(d_ref[i]));
  }
  Eigen::MatrixXd Kq(Eigen::MatrixXd::Zero(dimv, dimv)), 
                  Kv(Eigen::MatrixXd::Zero(dimv, dimv));
  for (int i=0; i<N; ++i) {
    const auto& lqr_policy_ref = riccati_recursion.getLQRPolicy();
    EXPECT_TRUE(lqr_policy[i].K.isApprox(lqr_policy_ref[i].K));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}