#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/unconstr/unconstr_ocp.hpp"
#include "idocp/riccati/unconstr_riccati_factorizer.hpp"
#include "idocp/riccati/unconstr_riccati_recursion.hpp"

#include "robot_factory.hpp"
#include "direction_factory.hpp"


namespace idocp {

class UnconstrRiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dimv = robot.dimv();
    dimx = 2 * robot.dimv();
    N = 20;
    T = 1;
    dt = T / N;

    kkt_matrix = KKTMatrix(robot, N);
    kkt_residual = KKTResidual(robot, N);
    d = Direction(robot, N);
    riccati_factorization = UnconstrRiccatiFactorization(N+1, SplitRiccatiFactorization(robot));
    for (int i=0; i<=N; ++i) {
      Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
      const Eigen::MatrixXd Q = seed * seed.transpose();
      kkt_matrix[i].Qxx = Q.topLeftCorner(dimx, dimx);
      kkt_matrix[i].Qxu = Q.topRightCorner(dimx, dimv);
      kkt_matrix[i].Qaa = Q.bottomRightCorner(dimv, dimv);
      kkt_residual[i].lx.setRandom();
      kkt_residual[i].la.setRandom();
      d[i].setRandom();
      seed = Eigen::MatrixXd::Random(dimv, dimv);
      riccati_factorization[i].Pqq = seed * seed.transpose();
      riccati_factorization[i].Pqv = Eigen::MatrixXd::Random(dimv, dimv);
      riccati_factorization[i].Pvq = riccati_factorization[i].Pqv.transpose();
      seed = Eigen::MatrixXd::Random(dimv, dimv);
      riccati_factorization[i].Pvv = seed * seed.transpose();
      riccati_factorization[i].sq.setRandom();
      riccati_factorization[i].sv.setRandom();
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
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  KKTResidual kkt_residual_ref = kkt_residual;
  UnconstrRiccatiRecursion riccati_recursion(robot, T, N);
  riccati_recursion.backwardRiccatiRecursion(kkt_matrix, kkt_residual, 
                                             riccati_factorization);
  riccati_factorization_ref[N].Pqq = kkt_matrix_ref[N].Qqq();
  riccati_factorization_ref[N].Pvv = kkt_matrix_ref[N].Qvv();
  riccati_factorization_ref[N].sq  = - kkt_residual_ref[N].lq();
  riccati_factorization_ref[N].sv  = - kkt_residual_ref[N].lv();
  EXPECT_TRUE(riccati_factorization[N].isApprox(riccati_factorization_ref[N]));
  UnconstrRiccatiFactorizer factorizer(robot);
  std::vector<LQRPolicy> lqr_policy(N, LQRPolicy(robot));
  for (int i=N-1; i>=0; --i) {
    factorizer.backwardRiccatiRecursion(
        riccati_factorization_ref[i+1], dt, 
        kkt_matrix_ref[i], kkt_residual_ref[i], riccati_factorization_ref[i],
        lqr_policy[i]);
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(kkt_matrix[i].isApprox(kkt_matrix_ref[i]));
    EXPECT_TRUE(kkt_residual[i].isApprox(kkt_residual_ref[i]));
    EXPECT_TRUE(riccati_factorization[i].isApprox(riccati_factorization_ref[i]));
  }
  d[0].dq().setRandom();
  d[0].dv().setRandom();
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(kkt_residual, d);
  for (int i=0; i<N; ++i) {
    factorizer.forwardRiccatiRecursion(kkt_residual_ref[i], dt, lqr_policy[i],
                                       d_ref[i], d_ref[i+1]);
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(d[i].isApprox(d_ref[i]));
  }
  Eigen::MatrixXd Kq(Eigen::MatrixXd::Zero(dimv, dimv)), 
                  Kv(Eigen::MatrixXd::Zero(dimv, dimv));
  for (int i=0; i<N; ++i) {
    riccati_recursion.getStateFeedbackGain(i, Kq, Kv);
    EXPECT_TRUE(Kq.isApprox(lqr_policy[i].Kq()));
    EXPECT_TRUE(Kv.isApprox(lqr_policy[i].Kv()));
  }

}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}