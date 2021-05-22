#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/parnmpc/unconstr_kkt_matrix_inverter.hpp"
#include "idocp/parnmpc/unconstr_split_backward_correction.hpp"

#include "robot_factory.hpp"


namespace idocp {

class UnconstrSplitBackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();
    dimx = 2*robot.dimv();
    dimKKT = 5*robot.dimv();

    kkt_matrix = SplitKKTMatrix(robot);
    const Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
    H = seed * seed.transpose();
    H.topRightCorner(dimv, dimx).setZero();
    H.block(dimv, dimx, dimv, dimv).setZero();
    kkt_matrix.Qaa = H.topLeftCorner(dimv, dimv);
    kkt_matrix.Qxu = H.bottomLeftCorner(dimx, dimv);
    kkt_matrix.Qxx = H.bottomRightCorner(dimx, dimx).transpose();

    kkt_residual = SplitKKTResidual(robot);
    kkt_residual.Fx.setRandom();
    kkt_residual.la.setRandom();
    kkt_residual.lx.setRandom();

    s = SplitSolution::Random(robot);
    d = SplitDirection::Random(robot);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
  int dimv, dimx, dimKKT;
  Eigen::MatrixXd H;
  SplitKKTMatrix kkt_matrix;
  SplitKKTResidual kkt_residual;
  SplitSolution s;
  SplitDirection d;
};


TEST_F(UnconstrSplitBackwardCorrectionTest, test) {
  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  UnconstrSplitBackwardCorrection corr(robot);
  auto s_new = SplitSolution::Random(robot);
  auto s_new_ref = s_new;
  corr.coarseUpdate(aux_mat_next, dt, kkt_matrix, kkt_residual, s, s_new);

  H.bottomRightCorner(dimx, dimx) += aux_mat_next;
  UnconstrKKTMatrixInverter inverter(robot);
  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  inverter.invert(dt, H, KKT_mat_inv);
  Eigen::VectorXd KKT_residual(Eigen::VectorXd::Zero(dimKKT));
  KKT_residual.head(dimx) = kkt_residual.Fx;
  KKT_residual.segment(dimx, dimv) = kkt_residual.la;
  KKT_residual.tail(dimx) = kkt_residual.lx;
  Eigen::VectorXd d_coarse = KKT_mat_inv * KKT_residual;
  s_new_ref.lmd = s.lmd - d_coarse.head(dimv);
  s_new_ref.gmm = s.gmm - d_coarse.segment(dimv, dimv);
  s_new_ref.a   = s.a   - d_coarse.segment(2*dimv, dimv);
  s_new_ref.q   = s.q   - d_coarse.segment(3*dimv, dimv);
  s_new_ref.v   = s.v   - d_coarse.segment(4*dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));

  const auto s_prev = SplitSolution::Random(robot);
  const auto s_new_prev = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  const auto s_new_next = SplitSolution::Random(robot);
  corr.backwardCorrectionSerial(s_next, s_new_next, s_new);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  x_res.head(dimv) = s_new_next.lmd - s_next.lmd;
  x_res.tail(dimv) = s_new_next.gmm - s_next.gmm;
  Eigen::VectorXd dx = KKT_mat_inv.block(0, dimKKT-dimx, dimx, dimx) * x_res;
  s_new_ref.lmd -= dx.head(dimv);
  s_new_ref.gmm -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.backwardCorrectionParallel(s_new);
  d_coarse.tail(dimKKT-dimx)
      = KKT_mat_inv.block(dimx, dimKKT-dimx, dimKKT-dimx, dimx) * x_res;
  s_new_ref.a -= d_coarse.segment(2*dimv, dimv);
  s_new_ref.q -= d_coarse.segment(3*dimv, dimv);
  s_new_ref.v -= d_coarse.segment(4*dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionSerial(s_prev, s_new_prev, s_new);
  x_res.head(dimv) = s_new_prev.q - s_prev.q;
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  dx = KKT_mat_inv.block(dimKKT-dimx, 0, dimx, dimx) * x_res;
  s_new_ref.q -= dx.head(dimv);
  s_new_ref.v -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionParallel(s_new);
  d_coarse.head(dimKKT-dimx).noalias()
      = KKT_mat_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_ref.lmd -= d_coarse.head(dimv);
  s_new_ref.gmm -= d_coarse.segment(dimv, dimv);
  s_new_ref.a   -= d_coarse.segment(2*dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  auto d_ref = d;
  corr.computeDirection(s, s_new, d);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.da() = s_new_ref.a - s.a;
  d_ref.dq() = s_new_ref.q - s.q;
  d_ref.dv() = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(UnconstrSplitBackwardCorrectionTest, testTerminal) {
  UnconstrSplitBackwardCorrection corr(robot);
  auto s_new = SplitSolution::Random(robot);
  auto s_new_ref = s_new;
  corr.coarseUpdate(dt, kkt_matrix, kkt_residual, s, s_new);

  UnconstrKKTMatrixInverter inverter(robot);
  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  inverter.invert(dt, H, KKT_mat_inv);
  Eigen::VectorXd KKT_residual(Eigen::VectorXd::Zero(dimKKT));
  KKT_residual.head(dimx) = kkt_residual.Fx;
  KKT_residual.segment(dimx, dimv) = kkt_residual.la;
  KKT_residual.tail(dimx) = kkt_residual.lx;
  Eigen::VectorXd d_coarse = KKT_mat_inv * KKT_residual;
  s_new_ref.lmd = s.lmd - d_coarse.head(dimv);
  s_new_ref.gmm = s.gmm - d_coarse.segment(dimv, dimv);
  s_new_ref.a   = s.a   - d_coarse.segment(2*dimv, dimv);
  s_new_ref.q   = s.q   - d_coarse.segment(3*dimv, dimv);
  s_new_ref.v   = s.v   - d_coarse.segment(4*dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}