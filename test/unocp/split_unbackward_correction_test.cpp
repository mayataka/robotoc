#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/unocp/split_unbackward_correction.hpp"


namespace idocp {

class SplitUnBackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    robot = Robot(urdf);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    dimv = robot.dimv();
    dimx = 2*robot.dimv();
    dimKKT = 5*robot.dimv();

    unkkt_matrix = SplitUnKKTMatrix(robot);
    const Eigen::MatrixXd seed = Eigen::MatrixXd::Random(3*dimv, 3*dimv);
    unkkt_matrix.Q = seed * seed.transpose();

    unkkt_residual = SplitUnKKTResidual(robot);
    unkkt_residual.KKT_residual.setRandom();

    s = SplitSolution::Random(robot);
    d = SplitDirection::Random(robot);
  }

  virtual void TearDown() {
  }

  std::string urdf;
  Robot robot;
  double dtau;
  int dimv, dimx, dimKKT;
  SplitUnKKTMatrix unkkt_matrix;
  SplitUnKKTResidual unkkt_residual;
  SplitSolution s;
  SplitDirection d;
};


TEST_F(SplitUnBackwardCorrectionTest, test) {
  SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  const SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  SplitUnBackwardCorrection corr(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  SplitDirection d_ref = d;
  corr.coarseUpdate(aux_mat_next, dtau, unkkt_matrix, unkkt_residual, s, d, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(5*dimv, 5*dimv));
  unkkt_matrix_ref.Qxx() += aux_mat_next;
  unkkt_matrix_ref.invert(dtau, KKT_mat_inv);
  d_ref.split_direction = KKT_mat_inv * unkkt_residual.KKT_residual;
  s_new_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_ref.a   = s.a - d_ref.da();
  s_new_ref.q   = s.q - d_ref.dq();
  s_new_ref.v   = s.v - d_ref.dv();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));

  const SplitSolution s_prev = SplitSolution::Random(robot);
  const SplitSolution s_new_prev = SplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  const SplitSolution s_new_next = SplitSolution::Random(robot);
  corr.backwardCorrectionSerial(s_next, s_new_next, s_new);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  x_res.head(dimv) = s_new_next.lmd - s_next.lmd;
  x_res.tail(dimv) = s_new_next.gmm - s_next.gmm;
  Eigen::VectorXd dx = KKT_mat_inv.block(0, dimKKT-dimx, dimx, dimx) * x_res;
  s_new_ref.lmd -= dx.head(dimv);
  s_new_ref.gmm -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.backwardCorrectionParallel(d, s_new);
  d_ref.split_direction.tail(dimKKT-dimx)
      = KKT_mat_inv.block(dimx, dimKKT-dimx, dimKKT-dimx, dimx) * x_res;
  s_new_ref.a -= d_ref.da();
  s_new_ref.q -= d_ref.dq();
  s_new_ref.v -= d_ref.dv();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));

  corr.forwardCorrectionSerial(s_prev, s_new_prev, s_new);
  x_res.head(dimv) = s_new_prev.q - s_prev.q;
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  dx = KKT_mat_inv.block(dimKKT-dimx, 0, dimx, dimx) * x_res;
  s_new_ref.q -= dx.head(dimv);
  s_new_ref.v -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionParallel(d, s_new);
  d_ref.split_direction.head(dimKKT-dimx).noalias()
      = KKT_mat_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_ref.lmd -= d_ref.dlmd();
  s_new_ref.gmm -= d_ref.dgmm();
  s_new_ref.a   -= d_ref.da();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));

  corr.computeDirection(s, s_new, d);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.da() = s_new_ref.a - s.a;
  d_ref.dq() = s_new_ref.q - s.q;
  d_ref.dv() = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(SplitUnBackwardCorrectionTest, testTerminal) {
  SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  const SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  SplitUnBackwardCorrection corr(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  SplitDirection d_ref = d;
  corr.coarseUpdate(dtau, unkkt_matrix, unkkt_residual, s, d, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(5*dimv, 5*dimv));
  unkkt_matrix_ref.invert(dtau, KKT_mat_inv);
  d_ref.split_direction = KKT_mat_inv * unkkt_residual.KKT_residual;
  s_new_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_ref.a   = s.a - d_ref.da();
  s_new_ref.q   = s.q - d_ref.dq();
  s_new_ref.v   = s.v - d_ref.dv();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}