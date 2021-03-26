#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_backward_correction.hpp"

#include "robot_factory.hpp"
#include "kkt_factory.hpp"


namespace idocp {

class SplitBackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((signed int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;
  void testWithSwitchingConstraint(const Robot& robot) const;
  void testTerminal(const Robot& robot) const;

  double dt;
};


void SplitBackwardCorrectionTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimKKT = 4*robot.dimv() + robot.dimu();

  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_residual_ref = kkt_residual;

  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  SplitBackwardCorrection corr(robot);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  std::cout << "before coarseUpdate()" << std::endl;
  auto kkt_mat_tmp = kkt_matrix;
  kkt_mat_tmp.Qvq() = kkt_mat_tmp.Qqv().transpose();
  kkt_mat_tmp.Qux() = kkt_mat_tmp.Qxu().transpose();
  std::cout << "kkt_mat_tmp.Qss()" << std::endl;
  std::cout << kkt_mat_tmp.Qss() << std::endl;
  Eigen::LLT<Eigen::MatrixXd> llt_ref(kkt_mat_tmp.Qss());
  ASSERT_TRUE(llt_ref.info() == Eigen::Success);

  corr.coarseUpdate(robot, dt, aux_mat_next, kkt_matrix, kkt_residual, s, s_new);
  std::cout << "after coarseUpdate()" << std::endl;

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qux() = kkt_matrix_ref.Qxu().transpose();
  kkt_matrix_ref.Qxx() += aux_mat_next;
  SplitKKTMatrixInverter inverter(robot);
  std::cout << "before invert()" << std::endl;
  inverter.invert(dt, kkt_matrix.Jac(), kkt_matrix.Qss(), KKT_mat_inv);
  std::cout << "after invert()" << std::endl;
  Eigen::VectorXd d0_ref = KKT_mat_inv * kkt_residual.splitKKTResidual();
  s_new_ref.lmd = s.lmd - d0_ref.head(dimv);
  s_new_ref.gmm = s.gmm - d0_ref.segment(dimv, dimv);
  s_new_ref.u   = s.u - d0_ref.segment(dimx, dimu);
  robot.integrateConfiguration(s.q, d0_ref.segment(dimx+dimu, dimv), -1, s_new_ref.q);
  s_new_ref.v   = s.v - d0_ref.segment(dimx+dimu+dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
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

  corr.backwardCorrectionParallel(robot, s_new);
  d0_ref.tail(dimKKT-dimx)
      = KKT_mat_inv.block(dimx, dimKKT-dimx, dimKKT-dimx, dimx) * x_res;
  s_new_ref.u -= d0_ref.segment(dimx, dimu);
  robot.integrateConfiguration(d0_ref.segment(dimx+dimu, dimv), -1, s_new_ref.q);
  s_new_ref.v -= d0_ref.segment(dimx+dimu+dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionSerial(robot, s_prev, s_new_prev, s_new);
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res.head(dimv));
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  dx = KKT_mat_inv.block(dimKKT-dimx, 0, dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(dimv), -1, s_new_ref.q);
  s_new_ref.v -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionParallel(s_new);
  d0_ref.head(dimKKT-dimx)
      = KKT_mat_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_ref.lmd -= d0_ref.head(dimv);
  s_new_ref.gmm -= d0_ref.segment(dimv, dimv);
  s_new_ref.u   -= d0_ref.segment(dimx, dimu);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  corr.computeDirection(robot, s, s_new, d);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.du() = s_new_ref.u - s.u;
  robot.subtractConfiguration(s_new_ref.q, s.q, d_ref.dq());
  d_ref.dv() = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void SplitBackwardCorrectionTest::testWithSwitchingConstraint(const Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  if (robot.hasFloatingBase()) {
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  else {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  const int dimQ = 2*robot.dimv() + robot.dimu();
  const int dimi = impulse_status.dimf();
  const int dimKKT = 4*robot.dimv() + robot.dimu() + dimi;
  ASSERT_TRUE(dimi > 0);
  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_matrix.Pq().setRandom();
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.P().setRandom();
  auto kkt_residual_ref = kkt_residual;

  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  SplitBackwardCorrection corr(robot);
  SplitSolution s = SplitSolution::Random(robot, impulse_status);
  SplitSolution s_new = SplitSolution::Random(robot, impulse_status);
  SplitSolution s_new_ref = s_new;
  const SplitSolution s_next = SplitSolution::Random(robot);
  const SplitSolution s_new_next = SplitSolution::Random(robot);
  corr.coarseUpdate(robot, dt, aux_mat_next, kkt_matrix, kkt_residual, s, s_new);
  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qux() = kkt_matrix_ref.Qxu().transpose();
  kkt_matrix_ref.Qxx() += aux_mat_next;
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(dt, kkt_matrix.Jac(), kkt_matrix.Pq(), kkt_matrix.Qss(), KKT_mat_inv);
  Eigen::VectorXd d0_ref = KKT_mat_inv * kkt_residual.splitKKTResidual();
  s_new_ref.lmd = s.lmd - d0_ref.head(dimv);
  s_new_ref.gmm = s.gmm - d0_ref.segment(dimv, dimv);
  s_new_ref.xi_stack() = s.xi_stack() - d0_ref.segment(dimx, dimi);
  s_new_ref.u   = s.u - d0_ref.segment(dimx+dimi, dimu);
  robot.integrateConfiguration(s.q, d0_ref.segment(dimx+dimi+dimu, dimv), -1, s_new_ref.q);
  s_new_ref.v   = s.v - d0_ref.segment(dimx+dimi+dimu+dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));

  const SplitSolution s_prev = SplitSolution::Random(robot);
  const SplitSolution s_new_prev = SplitSolution::Random(robot);
  corr.backwardCorrectionSerial(s_next, s_new_next, s_new);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  x_res.head(dimv) = s_new_next.lmd - s_next.lmd;
  x_res.tail(dimv) = s_new_next.gmm - s_next.gmm;
  Eigen::VectorXd dx = KKT_mat_inv.block(0, dimKKT-dimx, dimx, dimx) * x_res;
  s_new_ref.lmd -= dx.head(dimv);
  s_new_ref.gmm -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.backwardCorrectionParallel(robot, s_new);
  d0_ref.tail(dimKKT-dimx)
      = KKT_mat_inv.block(dimx, dimKKT-dimx, dimKKT-dimx, dimx) * x_res;
  s_new_ref.xi_stack() -= d0_ref.segment(dimx, dimi);
  s_new_ref.u -= d0_ref.segment(dimx+dimi, dimu);
  robot.integrateConfiguration(d0_ref.segment(dimx+dimi+dimu, dimv), -1, s_new_ref.q);
  s_new_ref.v -= d0_ref.segment(dimx+dimi+dimu+dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionSerial(robot, s_prev, s_new_prev, s_new);
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res.head(dimv));
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  dx = KKT_mat_inv.block(dimKKT-dimx, 0, dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(dimv), -1, s_new_ref.q);
  s_new_ref.v -= dx.tail(dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  corr.forwardCorrectionParallel(s_new);
  d0_ref.head(dimKKT-dimx).noalias()
      = KKT_mat_inv.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_ref.lmd -= d0_ref.head(dimv);
  s_new_ref.gmm -= d0_ref.segment(dimv, dimv);
  s_new_ref.xi_stack() -= d0_ref.segment(dimx, dimi);
  s_new_ref.u   -= d0_ref.segment(dimx+dimi, dimu);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  corr.computeDirection(robot, s, s_new, d);
  d_ref.setImpulseStatus(impulse_status);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.dxi() = s_new_ref.xi_stack() - s.xi_stack();
  d_ref.du() = s_new_ref.u - s.u;
  robot.subtractConfiguration(s_new_ref.q, s.q, d_ref.dq());
  d_ref.dv() = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void SplitBackwardCorrectionTest::testTerminal(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimKKT = 4*robot.dimv() + robot.dimu();

  auto kkt_matrix = testhelper::CreateSplitKKTMatrix(robot, dt);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = testhelper::CreateSplitKKTResidual(robot);
  auto kkt_residual_ref = kkt_residual;

  SplitBackwardCorrection corr(robot);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_new = SplitSolution::Random(robot);
  SplitSolution s_new_ref = s_new;
  corr.coarseUpdate(robot, dt, kkt_matrix, kkt_residual, s, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  kkt_matrix_ref.Qvq() = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qux() = kkt_matrix_ref.Qxu().transpose();
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(dt, kkt_matrix_ref.Jac(), kkt_matrix_ref.Qss(), KKT_mat_inv);
  Eigen::VectorXd d0_ref = KKT_mat_inv * kkt_residual.splitKKTResidual();
  s_new_ref.lmd = s.lmd - d0_ref.head(dimv);
  s_new_ref.gmm = s.gmm - d0_ref.segment(dimv, dimv);
  s_new_ref.u   = s.u - d0_ref.segment(dimx, dimu);
  robot.integrateConfiguration(s.q, d0_ref.segment(dimx+dimu, dimv), -1, s_new_ref.q);
  s_new_ref.v   = s.v - d0_ref.segment(dimx+dimu+dimv, dimv);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));
  EXPECT_TRUE(corr.auxMat().isApprox(KKT_mat_inv.topLeftCorner(dimx, dimx)));
}


TEST_F(SplitBackwardCorrectionTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  std::cout << "aaa" << std::endl;
  test(robot);
  std::cout << "bbb" << std::endl;
  testWithSwitchingConstraint(robot);
  std::cout << "ccc" << std::endl;
  testTerminal(robot);
  std::cout << "ddd" << std::endl;
}


TEST_F(SplitBackwardCorrectionTest, floating_base) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
  testWithSwitchingConstraint(robot);
  testTerminal(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}