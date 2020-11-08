#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/backward_correction.hpp"


namespace idocp {

class BackwardCorrectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testCoarseUpdate(const Robot& robot) const;
  void testBackwardCorrection(const Robot& robot) const;
  void testForkwardCorrection(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
  double dtau;
};


void BackwardCorrectionTest::testCoarseUpdate(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimKKT = 4*dimv + dimu;
  BackwardCorrection backward_correction(robot);
  const SplitSolution s = SplitSolution::Random(robot);
  SplitDirection d(robot);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxu().setRandom();
  kkt_matrix.Quu().setRandom();
  kkt_matrix.Fxx().setRandom();
  kkt_matrix.Fqq().setIdentity();
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvu().setRandom();
  kkt_residual.KKT_residual.setRandom();
  SplitSolution s_new_coarse(robot);
  backward_correction.coarseUpdate(robot, s, d, kkt_matrix, kkt_residual, s_new_coarse);
  Eigen::MatrixXd kkt_matrix_inverse = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_matrix.invert(kkt_matrix_inverse);
  SplitDirection d_ref(robot);
  d_ref.split_direction = kkt_matrix_inverse * kkt_residual.KKT_residual;
  SplitSolution s_new_coarse_ref(robot);
  s_new_coarse_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_coarse_ref.gmm = s.gmm - d_ref.dgmm();
  if (robot.has_floating_base()) {
    robot.integrateConfiguration(s.q, d_ref.dq(), -1, s_new_coarse_ref.q);
  }
  else {
    s_new_coarse_ref.q = s.q - d_ref.dq();
  }
  s_new_coarse_ref.v = s.v - d_ref.dv();
  s_new_coarse_ref.u = s.u - d_ref.du();
  EXPECT_TRUE(s_new_coarse.isApprox(s_new_coarse_ref));
  EXPECT_TRUE(d_ref.isApprox(d));
}


void BackwardCorrectionTest::testBackwardCorrection(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2 * robot.dimv();
  const int dimKKT = 4*dimv + dimu;
  BackwardCorrection backward_correction(robot);
  const SplitSolution s = SplitSolution::Random(robot);
  SplitDirection d(robot);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxu().setRandom();
  kkt_matrix.Quu().setRandom();
  kkt_matrix.Fxx().setRandom();
  kkt_matrix.Fqq().setIdentity();
  if (robot.has_floating_base()) {
    kkt_matrix.Fqq().topLeftCorner(6, 6).setRandom();
  }
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvu().setRandom();
  kkt_residual.KKT_residual.setRandom();
  SplitSolution s_new_coarse(robot);
  backward_correction.coarseUpdate(robot, s, d, kkt_matrix, kkt_residual, s_new_coarse);
  Eigen::MatrixXd kkt_matrix_inverse = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_matrix.invert(kkt_matrix_inverse);
  SplitDirection d_ref = d;
  SplitSolution s_new_coarse_ref = s_new_coarse;
  const SplitSolution s_next = SplitSolution::Random(robot);
  const SplitSolution s_new_next = SplitSolution::Random(robot);
  backward_correction.backwardCorrectionSerial(s_next, s_new_next, s_new_coarse);
  Eigen::VectorXd lmdgmm_res = Eigen::VectorXd::Zero(dimx);
  lmdgmm_res.head(dimv) = s_new_next.lmd - s_next.lmd;
  lmdgmm_res.tail(dimv) = s_new_next.gmm - s_next.gmm;
  const Eigen::VectorXd dlmdgmm = kkt_matrix_inverse.topRightCorner(dimx, dimx) * lmdgmm_res;
  s_new_coarse_ref.lmd -= dlmdgmm.head(dimv);
  s_new_coarse_ref.gmm -= dlmdgmm.tail(dimv);
  EXPECT_TRUE(s_new_coarse.isApprox(s_new_coarse_ref));
  backward_correction.backwardCorrectionParallel(robot, d, s_new_coarse);
  d_ref.split_direction.tail(dimKKT-dimx).noalias()
      = kkt_matrix_inverse.bottomRightCorner(dimKKT-dimx, dimx) * lmdgmm_res;
  s_new_coarse_ref.u.noalias() -= d_ref.du();
  robot.integrateConfiguration(d_ref.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= d_ref.dv();
  EXPECT_TRUE(s_new_coarse.isApprox(s_new_coarse_ref));
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const SplitSolution s_new_prev = SplitSolution::Random(robot);
  backward_correction.forwardCorrectionSerial(robot, s_prev, s_new_prev, s_new_coarse);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res.head(dimv));
  x_res.tail(dimv) = s_new_prev.v - s_prev.v;
  Eigen::VectorXd dx = Eigen::VectorXd::Zero(dimx);
  dx.noalias() = kkt_matrix_inverse.bottomLeftCorner(dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(dimv), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= dx.tail(dimv);
  EXPECT_TRUE(s_new_coarse.isApprox(s_new_coarse_ref));
  backward_correction.forwardCorrectionParallel(d_ref, s_new_coarse);
  d_ref.split_direction.head(dimKKT-dimx).noalias()
      = kkt_matrix_inverse.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_coarse_ref.lmd.noalias() -= d_ref.dlmd();
  s_new_coarse_ref.gmm.noalias() -= d_ref.dgmm();
  s_new_coarse_ref.u.noalias() -= d_ref.du();
  EXPECT_TRUE(s_new_coarse.isApprox(s_new_coarse_ref));
}


TEST_F(BackwardCorrectionTest, fixed_base) {
  testCoarseUpdate(fixed_base_robot);
  testBackwardCorrection(fixed_base_robot);
}


TEST_F(BackwardCorrectionTest, floating_base) {
  testCoarseUpdate(floating_base_robot);
  testBackwardCorrection(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}