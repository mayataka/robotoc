#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_dynamics_data.hpp"
#include "robotoc/impulse/impulse_dynamics.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ImpulseDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test_computeResidual(Robot& robot, const ImpulseStatus& impulse_status);
  static void test_linearize(Robot& robot, const ImpulseStatus& impulse_status);
  static void test_condense(Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseDynamicsTest::test_computeResidual(Robot& robot, const ImpulseStatus& impulse_status) {
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamics id(robot);
  id.evalImpulseDynamics(robot, impulse_status, s);
  const double l1norm = id.constraintViolation();
  const double squarednorm = id.KKTError();
  ImpulseDynamicsData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  double l1norm_ref = data.ImDC().lpNorm<1>();
  double squarednorm_ref = data.ImDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


void ImpulseDynamicsTest::test_linearize(Robot& robot, const ImpulseStatus& impulse_status) {
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamics id(robot);
  auto kkt_residual = ImpulseSplitKKTResidual::Random(robot, impulse_status);
  auto kkt_residual_ref = kkt_residual;
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_residual);
  const double l1norm = id.constraintViolation();
  const double squarednorm = id.KKTError();
  ImpulseDynamicsData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq(), data.dImDddv);
  robot.computeImpulseVelocityDerivatives(impulse_status, data.dCdq(), data.dCdv());
  kkt_residual_ref.lq().noalias() += data.dImDdq().transpose() * s.beta;
  kkt_residual_ref.ldv.noalias()  += data.dImDddv.transpose() * s.beta;
  kkt_residual_ref.lf().noalias() -= data.dCdv() * s.beta;
  kkt_residual_ref.lq().noalias() += data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv().noalias() += data.dCdv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv.noalias()  += data.dCdv().transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  double l1norm_ref = data.ImDC().lpNorm<1>();
  double squarednorm_ref = data.ImDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


void ImpulseDynamicsTest::test_condense(Robot& robot, const ImpulseStatus& impulse_status) {
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamics id(robot);
  auto kkt_residual = ImpulseSplitKKTResidual::Random(robot, impulse_status);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_residual);
  ImpulseDynamicsData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data.dImDdq(), data.dImDddv);
  robot.computeImpulseVelocityDerivatives(impulse_status, data.dCdq(), data.dCdv());
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimi();
  auto kkt_matrix = ImpulseSplitKKTMatrix::Random(robot, impulse_status);
  kkt_matrix.Fxx.setZero();
  kkt_matrix.Qdvdv.setZero();
  kkt_matrix.Qdvdv.diagonal().setRandom();
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix, kkt_residual);
  robot.computeMJtJinv(data.dImDddv, data.dCdv(), data.MJtJinv());
  data.MJtJinv_dImDCdqv() = data.MJtJinv() * data.dImDCdqv();
  data.MJtJinv_ImDC()     = data.MJtJinv() * data.ImDC();
  Eigen::MatrixXd Qdvdvff = Eigen::MatrixXd::Zero(dimv+dimf, dimv+dimf);
  Qdvdvff.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Qdvdv;
  Qdvdvff.bottomRightCorner(dimf, dimf) = kkt_matrix_ref.Qff();
  data.Qdvfqv() = - Qdvdvff * data.MJtJinv_dImDCdqv();
  data.Qdvfqv().bottomLeftCorner(dimf, dimv) -= kkt_matrix_ref.Qqf().transpose();
  data.ldvf().head(dimv) = kkt_residual_ref.ldv;
  data.ldvf().tail(dimf) = - kkt_residual_ref.lf();
  data.ldvf() -= Qdvdvff * data.MJtJinv_ImDC();
  kkt_matrix_ref.Qxx -= data.MJtJinv_dImDCdqv().transpose() * data.Qdvfqv();
  kkt_matrix_ref.Qxx.topRows(dimv) += kkt_matrix_ref.Qqf() * data.MJtJinv_dImDCdqv().bottomRows(dimf);
  kkt_residual_ref.lx -= data.MJtJinv_dImDCdqv().transpose() * data.ldvf();
  kkt_residual_ref.lx.head(dimv) += kkt_matrix_ref.Qqf() * data.MJtJinv_ImDC().tail(dimf);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv).setIdentity();
  kkt_matrix_ref.Fvv().setIdentity();
  kkt_matrix_ref.Fxx -= OOIO_mat * data.MJtJinv_dImDCdqv();
  kkt_residual_ref.Fx -= OOIO_mat * data.MJtJinv_ImDC();
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  auto d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  id.expandPrimal(d);
  d_ref.ddvf() = - data.MJtJinv() * (data.dImDCdqv() * d_ref.dx + data.ImDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));
  const auto d_next = ImpulseSplitDirection::Random(robot);
  id.expandDual(d_next, d);
  d_ref.dbetamu() = - data.MJtJinv() * (data.Qdvfqv() * d.dx 
                                         + OOIO_mat.transpose() * d_next.dlmdgmm
                                         + data.ldvf());
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ImpulseDynamicsTest, robotManipulator) {
  const double dt = 0.001;
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
}


TEST_F(ImpulseDynamicsTest, quadrupedalRobot) {
  const double dt = 0.001;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
}


TEST_F(ImpulseDynamicsTest, humanoidRobot) {
  const double dt = 0.001;
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test_computeResidual(robot, impulse_status);
  test_linearize(robot, impulse_status);
  test_condense(robot, impulse_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}