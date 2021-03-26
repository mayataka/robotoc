#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_backward_correction.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseSplitBackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((signed int) time(0));
  }

  virtual void TearDown() {
  }

  ImpulseSplitKKTMatrix createKKTMatrix(const Robot& robot, const ImpulseStatus& impulse_status) const;
  ImpulseSplitKKTResidual createKKTResidual(const Robot& robot, const ImpulseStatus& impulse_status) const;

  void test(const Robot& robot) const;
};


ImpulseSplitKKTMatrix ImpulseSplitBackwardCorrectionTest::createKKTMatrix(
    const Robot& robot, const ImpulseStatus& impulse_status) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimi = impulse_status.dimf();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(2*dimv+dimi, 2*dimv+dimi);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_matrix.Qss() = seed * seed.transpose();
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationMinus(robot.generateFeasibleConfiguration(), 
                                       robot.generateFeasibleConfiguration(), kkt_matrix.Fqq());
  }
  kkt_matrix.Fvf().setRandom();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Vq().setRandom();
  kkt_matrix.Vv().setRandom();
  return kkt_matrix;
}


ImpulseSplitKKTResidual ImpulseSplitBackwardCorrectionTest::createKKTResidual(
    const Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.splitKKTResidual().setRandom();
  return kkt_residual;
}


void ImpulseSplitBackwardCorrectionTest::test(const Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimi = impulse_status.dimf();
  const int dimKKT = 4*robot.dimv() + 2*dimi;

  auto kkt_matrix = createKKTMatrix(robot, impulse_status);
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual = createKKTResidual(robot, impulse_status);
  auto kkt_residual_ref = kkt_residual;

  Eigen::MatrixXd aux_mat_seed(Eigen::MatrixXd::Random(dimx, dimx));
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed * aux_mat_seed.transpose();
  ImpulseSplitBackwardCorrection corr(robot);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitSolution s_new = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitSolution s_new_ref = s_new;
  corr.coarseUpdate(robot, aux_mat_next, kkt_matrix, kkt_residual, s, s_new);

  Eigen::MatrixXd KKT_mat_inv(Eigen::MatrixXd::Zero(dimKKT, dimKKT));
  ImpulseSplitKKTMatrixInverter inverter(robot);
  kkt_matrix_ref.Qxx() += aux_mat_next;
  inverter.invert(kkt_matrix_ref.Jac(), kkt_matrix_ref.Qss(), KKT_mat_inv);
  Eigen::VectorXd d0_ref = KKT_mat_inv * kkt_residual.splitKKTResidual();
  s_new_ref.lmd = s.lmd - d0_ref.head(dimv);
  s_new_ref.gmm = s.gmm - d0_ref.segment(dimv, dimv);
  s_new_ref.mu_stack() = s.mu_stack() - d0_ref.segment(dimx, dimi);
  s_new_ref.f_stack()  = s.f_stack()  - d0_ref.segment(dimx+dimi, dimi);
  robot.integrateConfiguration(s.q, d0_ref.segment(dimx+2*dimi, dimv), -1, s_new_ref.q);
  s_new_ref.v   = s.v - d0_ref.segment(dimx+2*dimi+dimv, dimv);
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
  s_new_ref.mu_stack() -= d0_ref.segment(dimx, dimi);
  s_new_ref.f_stack()  -= d0_ref.segment(dimx+dimi, dimi);
  robot.integrateConfiguration(d0_ref.segment(dimx+2*dimi, dimv), -1, s_new_ref.q);
  s_new_ref.v -= d0_ref.segment(dimx+2*dimi+dimv, dimv);
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
  s_new_ref.mu_stack()-= d0_ref.segment(dimx, dimi);
  s_new_ref.f_stack() -= d0_ref.segment(dimx+dimi, dimi);
  EXPECT_TRUE(s_new.isApprox(s_new_ref));

  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitDirection d_ref = d;
  corr.computeDirection(robot, s, s_new, d);
  d_ref.dlmd() = s_new_ref.lmd - s.lmd;
  d_ref.dgmm() = s_new_ref.gmm - s.gmm;
  d_ref.dmu()  = s_new_ref.mu_stack() - s.mu_stack();
  d_ref.df()   = s_new_ref.f_stack() - s.f_stack();
  robot.subtractConfiguration(s_new_ref.q, s.q, d_ref.dq());
  d_ref.dv()   = s_new_ref.v - s.v;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ImpulseSplitBackwardCorrectionTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(ImpulseSplitBackwardCorrectionTest, floating_base) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}