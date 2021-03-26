#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler_data.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseDynamicsForwardEulerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void testLinearizeInverseImpulseDynamics(Robot& robot, 
                                                  const ImpulseStatus& impulse_status);
  static void testLinearizeImpulseVelocityConstraints(Robot& robot, 
                                                      const ImpulseStatus& impulse_status);
  static void testLinearizeImpulseDynamics(Robot& robot, 
                                           const ImpulseStatus& impulse_status);
  static void testCondensing(Robot& robot, const ImpulseStatus& impulse_status);
  static void testExpansionPrimal(Robot& robot, const ImpulseStatus& impulse_status);
  static void testExpansionDual(Robot& robot, const ImpulseStatus& impulse_status);
  static void testIntegration(Robot& robot, const ImpulseStatus& impulse_status);
  static void testComputeResidual(Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseDynamicsForwardEulerTest::testLinearizeInverseImpulseDynamics(Robot& robot, 
                                                                          const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsForwardEulerData data(robot), data_ref(robot);
  ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_ref.ImD());
  robot.RNEAImpulseDerivatives(s.q, s.dv, data_ref.dImDdq(), data_ref.dImDddv);
  EXPECT_TRUE(data_ref.ImDC().isApprox(data.ImDC()));
  EXPECT_TRUE(data_ref.dImDddv.isApprox(data.dImDddv));
  EXPECT_TRUE(data_ref.dImDdq().isApprox(data.dImDdq()));
}


void ImpulseDynamicsForwardEulerTest::testLinearizeImpulseVelocityConstraints(Robot& robot, 
                                                                              const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsForwardEulerData data(robot), data_ref(robot);
  data.setImpulseStatus(impulse_status);
  data_ref.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamicsForwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, data);
  robot.computeImpulseVelocityResidual(impulse_status, data_ref.C());
  robot.computeImpulseVelocityDerivatives(impulse_status, data_ref.dCdq(), 
                                          data_ref.dCdv());
  EXPECT_TRUE(data.ImDC().isApprox(data_ref.ImDC()));
  EXPECT_TRUE(data.dCdq().isApprox(data_ref.dCdq()));
  EXPECT_TRUE(data.dCdv().isApprox(data_ref.dCdv()));
  EXPECT_TRUE(data.dImDCdqv().isApprox(data_ref.dImDCdqv()));
}


void ImpulseDynamicsForwardEulerTest::testLinearizeImpulseDynamics(Robot& robot, 
                                                                   const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom();
  kkt_residual.ldv.setRandom();
  kkt_residual.lf().setRandom();
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, s.v+s.dv);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data);
  ImpulseDynamicsForwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, data);
  Eigen::MatrixXd dImDdf = Eigen::MatrixXd::Zero(robot.dimv(), impulse_status.dimf());
  robot.updateKinematics(s.q, s.v+s.dv);
  auto contact_status = robot.createContactStatus();
  contact_status.setActivity(impulse_status.isImpulseActive());
  robot.dRNEAPartialdFext(contact_status, dImDdf);
  kkt_residual_ref.lq() += data.dImDdq().transpose() * s.beta + data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += data.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.ldv += data.dImDddv.transpose() * s.beta + data.dCdv().transpose() * s.mu_stack();
  kkt_residual_ref.lf() += dImDdf.transpose() * s.beta;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  const double l1norm_ref = data.ImDC().lpNorm<1>();
  const double squarednorm_ref = data.ImDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm_ref, l1norm);
  EXPECT_DOUBLE_EQ(squarednorm_ref, squarednorm);
}  


void ImpulseDynamicsForwardEulerTest::testCondensing(Robot& robot, 
                                                     const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.lx().setRandom();
  kkt_residual.ldv.setRandom();
  kkt_residual.lf().setRandom();
  kkt_residual.Fx().setRandom();
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qdvdv().diagonal().setRandom();
  const Eigen::MatrixXd Qff_seed = Eigen::MatrixXd::Random(dimf, dimf);
  kkt_matrix.Qff() = Qff_seed * Qff_seed.transpose();
  kkt_matrix.Qqf().setRandom();
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseDynamicsForwardEuler id(robot);
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  data.dImDddv.setRandom();
  data.dImDCdqv().setRandom();
  data.dImDCdqv().topRightCorner(dimv, dimv).setZero();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dImDCdqv().setRandom();
  data.Qdvfqv().setRandom();
  data.ImDC().setRandom();
  data.MJtJinv_ImDC().setRandom();
  data.ldvf().setRandom();
  ImpulseDynamicsForwardEulerData data_ref = data;
  data_ref.MJtJinv_dImDCdqv().setZero();
  data_ref.Qdvfqv().setZero();
  data_ref.MJtJinv_ImDC().setZero();
  data_ref.ldvf().setZero();
  ImpulseDynamicsForwardEuler::condensing(robot, impulse_status, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dImDCdqv() = data_ref.MJtJinv() * data_ref.dImDCdqv();
  data_ref.MJtJinv_ImDC() = data_ref.MJtJinv() * data_ref.ImDC();
  Eigen::MatrixXd Qdvdvff = Eigen::MatrixXd::Zero(dimv+dimf, dimv+dimf);
  Qdvdvff.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Qdvdv();
  Qdvdvff.bottomRightCorner(dimf, dimf) = kkt_matrix_ref.Qff();
  data_ref.Qdvfqv() = - Qdvdvff * data_ref.MJtJinv_dImDCdqv();
  data_ref.Qdvfqv().bottomLeftCorner(dimf, dimv) -= kkt_matrix_ref.Qqf().transpose();
  data_ref.ldvf().head(dimv) = kkt_residual_ref.ldv;
  data_ref.ldvf().tail(dimf) = - kkt_residual_ref.lf();
  data_ref.ldvf() -= Qdvdvff * data_ref.MJtJinv() * data_ref.ImDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dImDCdqv().transpose() * data_ref.Qdvfqv();
  kkt_matrix_ref.Qxx().topRows(dimv) += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_dImDCdqv().bottomRows(dimf);
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dImDCdqv().transpose() * data_ref.ldvf();
  kkt_residual_ref.lx().head(dimv) += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_ImDC().tail(dimf);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv).setIdentity();
  kkt_matrix_ref.Fvv().setIdentity();
  kkt_matrix_ref.Fxx() -= OOIO_mat * data_ref.MJtJinv_dImDCdqv();
  kkt_residual_ref.Fx() -= OOIO_mat * data_ref.MJtJinv_ImDC();
  EXPECT_TRUE(data_ref.MJtJinv().isApprox(data.MJtJinv()));
  EXPECT_TRUE(data_ref.MJtJinv_dImDCdqv().isApprox(data.MJtJinv_dImDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_ImDC().isApprox(data.MJtJinv_ImDC()));
  EXPECT_TRUE(data_ref.Qdvfqv().isApprox(data.Qdvfqv()));
  EXPECT_TRUE(data_ref.ldvf().isApprox(data.ldvf()));
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
}


void ImpulseDynamicsForwardEulerTest::testExpansionPrimal(Robot& robot, 
                                                          const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  data.dImDddv.setRandom();
  data.dImDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.ImDC().setRandom();
  data.MJtJinv_dImDCdqv() = data.MJtJinv() * data.dImDCdqv();
  data.MJtJinv_ImDC() = data.MJtJinv() * data.ImDC();
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitDirection d_ref = d;
  ImpulseDynamicsForwardEuler::expansionPrimal(robot, data, d);
  d_ref.ddvf() = - data.MJtJinv() * (data.dImDCdqv() * d.dx() + data.ImDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ImpulseDynamicsForwardEulerTest::testExpansionDual(Robot& robot, 
                                                        const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  data.dImDddv.setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.Qdvfqv().setRandom();
  data.ldvf().setRandom();
  const ImpulseDynamicsForwardEulerData data_ref = data;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*dimv);
  const Eigen::VectorXd dgmm_next = dlmdgmm_next.tail(dimv);
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitDirection d_ref = d;
  ImpulseDynamicsForwardEuler::expansionDual(robot, data, kkt_matrix, kkt_residual, dgmm_next, d);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv).setIdentity();
  d_ref.dbetamu() = - data_ref.MJtJinv() * (data_ref.Qdvfqv() * d.dx() 
                                            + OOIO_mat.transpose() * dlmdgmm_next 
                                            + data_ref.ldvf());
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ImpulseDynamicsForwardEulerTest::testIntegration(Robot& robot, 
                                                      const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.lx().setRandom();
  kkt_residual.ldv.setRandom();
  kkt_residual.Fx().setRandom();
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qdvdv().diagonal().setRandom();
  const Eigen::MatrixXd Qff_seed = Eigen::MatrixXd::Random(dimf, dimf);
  kkt_matrix.Qff() = Qff_seed * Qff_seed.transpose();
  kkt_matrix.Qqf().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  robot.updateKinematics(s.q, s.v+s.dv);
  ImpulseDynamicsForwardEuler id(robot), id_ref(robot);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix, kkt_residual);
  ImpulseDynamicsForwardEulerData data_ref(robot);
  data_ref.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  id_ref.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsForwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data_ref);
  ImpulseDynamicsForwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, data_ref);
  robot.computeMJtJinv(data_ref.dImDddv, data_ref.dCdv(), data_ref.MJtJinv());
  ImpulseDynamicsForwardEuler::condensing(robot, impulse_status, data_ref, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitDirection d_ref = d;
  id.computeCondensedPrimalDirection(robot, d);
  ImpulseDynamicsForwardEuler::expansionPrimal(robot, data_ref, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  const Eigen::VectorXd dgmm_next = Eigen::VectorXd::Random(dimv);
  id.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, dgmm_next, d);
  ImpulseDynamicsForwardEuler::expansionDual(robot, data_ref, kkt_matrix_ref, kkt_residual_ref, dgmm_next, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ImpulseDynamicsForwardEulerTest::testComputeResidual(Robot& robot, 
                                                          const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsForwardEuler id(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  robot.updateKinematics(s.q, s.v+s.dv);
  id.computeImpulseDynamicsResidual(robot, impulse_status, s, kkt_residual);
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD());
  robot.updateKinematics(s.q, s.v+s.dv);
  robot.computeImpulseVelocityResidual(impulse_status, data.C());
  double l1norm_ref = data.ImDC().lpNorm<1>();
  double squarednorm_ref = data.ImDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ImpulseDynamicsForwardEulerTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testExpansionPrimal(robot, impulse_status);
  testExpansionDual(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
  impulse_status.activateImpulse(0);
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testExpansionPrimal(robot, impulse_status);
  testExpansionDual(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
}


TEST_F(ImpulseDynamicsForwardEulerTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testExpansionPrimal(robot, impulse_status);
  testExpansionDual(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testExpansionPrimal(robot, impulse_status);
  testExpansionDual(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}