#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_backward_correction.hpp"

#include "test_helper.hpp"

namespace idocp {

class SplitBackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((signed int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  SplitKKTMatrix createKKTMatrix(const Robot& robot) const;
  SplitKKTResidual createKKTResidual(const Robot& robot) const;

  void test(const Robot& robot) const;
  void testTerminal(const Robot& robot) const;

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


SplitKKTMatrix SplitBackwardCorrectionTest::createKKTMatrix(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(2*dimv+dimu, 2*dimv+dimu);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qss() = seed * seed.transpose();
  kkt_matrix.Qvq().setZero();
  kkt_matrix.Quq().setZero();
  kkt_matrix.Quv().setZero();
  kkt_matrix.Fqq() = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  return kkt_matrix;
}


SplitKKTResidual SplitBackwardCorrectionTest::createKKTResidual(const Robot& robot) const {
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


void SplitBackwardCorrectionTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimKKT = 4*robot.dimv() + robot.dimu();

  auto kkt_matrix = createKKTMatrix(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = createKKTResidual(robot);
  auto kkt_residual_ref = kkt_residual;

  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  SplitBackwardCorrection corr(robot);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  corr.coarseUpdate(robot, aux_mat_next, kkt_matrix, kkt_residual, s, d, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qux() = kkt_matrix_ref.Qxu().transpose();
  kkt_matrix_ref.Qxx() += aux_mat_next;
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(kkt_matrix.Jac(), kkt_matrix.Qss(), KKT_mat_inv);
  d_ref.split_direction = KKT_mat_inv * kkt_residual.KKT_residual;
  s_new_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_ref.u   = s.u - d_ref.du();
  robot.integrateConfiguration(s.q, d_ref.dq(), -1, s_new_ref.q);
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

  corr.backwardCorrectionParallel(robot, d, s_new);
  d_ref.split_direction.tail(dimKKT-dimx)
      = KKT_mat_inv.block(dimx, dimKKT-dimx, dimKKT-dimx, dimx) * x_res;
  s_new_ref.u -= d_ref.du();
  robot.integrateConfiguration(d_ref.dq(), -1, s_new_ref.q);
  s_new_ref.v -= d_ref.dv();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));

  corr.forwardCorrectionSerial(robot, s_prev, s_new_prev, s_new);
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res.head(dimv));
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  dx = KKT_mat_inv.block(dimKKT-dimx, 0, dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(dimv), -1, s_new_ref.q);
  s_new_ref.v -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionParallel(d, s_new);
  d_ref.split_direction.head(dimKKT-dimx).noalias()
      = KKT_mat_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_ref.lmd -= d_ref.dlmd();
  s_new_ref.gmm -= d_ref.dgmm();
  s_new_ref.u   -= d_ref.du();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));

  corr.computeDirection(robot, s, s_new, d);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.du() = s_new_ref.u - s.u;
  robot.subtractConfiguration(s_new_ref.q, s.q, d_ref.dq());
  d_ref.dv() = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void SplitBackwardCorrectionTest::testTerminal(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimKKT = 4*robot.dimv() + robot.dimu();

  auto kkt_matrix = createKKTMatrix(robot);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = createKKTResidual(robot);
  auto kkt_residual_ref = kkt_residual;
  SplitBackwardCorrection corr(robot);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  corr.coarseUpdate(robot, kkt_matrix, kkt_residual, s, d, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qux() = kkt_matrix_ref.Qxu().transpose();
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(kkt_matrix_ref.Jac(), kkt_matrix_ref.Qss(), KKT_mat_inv);
  d_ref.split_direction = KKT_mat_inv * kkt_residual.KKT_residual;
  s_new_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_ref.u   = s.u - d_ref.du();
  robot.integrateConfiguration(s.q, d_ref.dq(), -1, s_new_ref.q);
  s_new_ref.v   = s.v - d_ref.dv();
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));
}


TEST_F(SplitBackwardCorrectionTest, fixedBase) {
  test(fixed_base_robot);
  testTerminal(fixed_base_robot);
}


TEST_F(SplitBackwardCorrectionTest, floating_base) {
  test(floating_base_robot);
  // testTerminal(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}