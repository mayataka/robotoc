#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"
#include "idocp/unocp/unconstrained_container.hpp"
#include "idocp/unocp/unriccati_recursion.hpp"

#include "robot_factory.hpp"


namespace idocp {

class UnRiccatiRecursionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dimv = robot.dimv();
    N = 20;
    T = 1;
    dt = T / N;

    unkkt_matrix = UnKKTMatrix(N+1, robot);
    unkkt_residual = UnKKTResidual(N+1, robot);
    d = UnDirection(N+1, robot);
    riccati_factorization = UnRiccatiFactorization(N+1, robot);
    for (int i=0; i<=N; ++i) {
      Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
      unkkt_matrix[i].Q = seed * seed.transpose();
      unkkt_residual[i].KKT_residual.setRandom();
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
    Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
    terminal_kkt_matrix = SplitKKTMatrix(robot);
    terminal_kkt_matrix.Qqq() = seed * seed.transpose();
    seed = Eigen::MatrixXd::Random(dimv, dimv);
    terminal_kkt_matrix.Qvv() = seed * seed.transpose();
    terminal_kkt_residual = SplitKKTResidual(robot);   
    terminal_kkt_residual.Fx.setRandom();
    terminal_kkt_residual.lx.setRandom();
  }

  virtual void TearDown() {
  }

  Robot robot;
  int N, dimv;
  double T, dt;
  UnKKTMatrix unkkt_matrix;
  UnKKTResidual unkkt_residual;
  SplitKKTMatrix terminal_kkt_matrix;
  SplitKKTResidual terminal_kkt_residual;
  UnDirection d;
  UnRiccatiFactorization riccati_factorization;
};


TEST_F(UnRiccatiRecursionTest, test) {
  UnRiccatiFactorization riccati_factorization_ref = riccati_factorization;
  UnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  UnKKTResidual unkkt_residual_ref = unkkt_residual;
  SplitKKTMatrix terminal_kkt_matrix_ref = terminal_kkt_matrix;
  SplitKKTResidual terminal_kkt_residual_ref = terminal_kkt_residual;
  UnRiccatiRecursion riccati_recursion(robot, T, N);
  riccati_recursion.backwardRiccatiRecursionTerminal(terminal_kkt_matrix, 
                                                     terminal_kkt_residual, 
                                                     riccati_factorization);
  riccati_factorization_ref[N].Pqq = terminal_kkt_matrix_ref.Qqq();
  riccati_factorization_ref[N].Pvv = terminal_kkt_matrix_ref.Qvv();
  riccati_factorization_ref[N].sq  = - terminal_kkt_residual_ref.lq();
  riccati_factorization_ref[N].sv  = - terminal_kkt_residual_ref.lv();
  EXPECT_TRUE(riccati_factorization[N].isApprox(riccati_factorization_ref[N]));
  riccati_recursion.backwardRiccatiRecursion(unkkt_matrix, unkkt_residual, 
                                             riccati_factorization);
  UnRiccatiFactorizer factorizer(N, SplitUnRiccatiFactorizer(robot));
  for (int i=N-1; i>=0; --i) {
    factorizer[i].backwardRiccatiRecursion(
        riccati_factorization_ref[i+1], dt, 
        unkkt_matrix_ref[i], unkkt_residual_ref[i], riccati_factorization_ref[i]);
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(unkkt_matrix[i].isApprox(unkkt_matrix_ref[i]));
    EXPECT_TRUE(unkkt_residual[i].isApprox(unkkt_residual_ref[i]));
    EXPECT_TRUE(riccati_factorization[i].isApprox(riccati_factorization_ref[i]));
  }
  d[0].dq().setRandom();
  d[0].dv().setRandom();
  UnDirection d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(unkkt_residual, d);
  for (int i=0; i<N; ++i) {
    factorizer[i].forwardRiccatiRecursion(unkkt_residual_ref[i], d_ref[i],  
                                          dt, d_ref[i+1]);
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(d[i].isApprox(d_ref[i]));
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}