#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

class ImpulseDynamicsBackwardEulerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void testLinearizeInverseImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status);
  static void testLinearizeImpulseVelocityConstraints(Robot& robot, const ImpulseStatus& impulse_status);
  static void testLinearizeImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status);
  static void testCondensing(Robot& robot, const ImpulseStatus& impulse_status);
  static void testIntegration(Robot& robot, const ImpulseStatus& impulse_status);
  static void testComputeResidual(Robot& robot, const ImpulseStatus& impulse_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseDynamicsBackwardEulerTest::testLinearizeInverseImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsBackwardEulerData data(robot), data_ref(robot);
  ImpulseDynamicsBackwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_ref.ImD);
  robot.RNEAImpulseDerivatives(s.q, s.dv, data_ref.dImDdq, data_ref.dImDddv);
  EXPECT_TRUE(data_ref.ImD.isApprox(data.ImD));
  EXPECT_TRUE(data_ref.dImDddv.isApprox(data.dImDddv));
  EXPECT_TRUE(data_ref.dImDdq.isApprox(data.dImDdq));
}


void ImpulseDynamicsBackwardEulerTest::testLinearizeImpulseVelocityConstraints(Robot& robot, const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseSplitKKTMatrix kkt_matrix(robot), kkt_matrix_ref(robot);
  ImpulseSplitKKTResidual kkt_residual(robot), kkt_residual_ref(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  ImpulseDynamicsBackwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, kkt_matrix, kkt_residual);
  robot.computeImpulseVelocityResidual(impulse_status, kkt_residual_ref.V());
  robot.computeImpulseVelocityDerivatives(impulse_status, kkt_matrix_ref.Vq(), kkt_matrix_ref.Vv());
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void ImpulseDynamicsBackwardEulerTest::testLinearizeImpulseDynamics(Robot& robot, const ImpulseStatus& impulse_status) {
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
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  ImpulseDynamicsBackwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  ImpulseDynamicsBackwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data);
  ImpulseDynamicsBackwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, kkt_matrix_ref, kkt_residual_ref);
  Eigen::MatrixXd dImDdf = Eigen::MatrixXd::Zero(robot.dimv(), impulse_status.dimf());
  robot.updateKinematics(s.q, s.v);
  ContactStatus contact_status(impulse_status.maxPointContacts());
  contact_status.setContactStatus(impulse_status.isImpulseActive());
  robot.dRNEAPartialdFext(contact_status, dImDdf);
  kkt_residual_ref.lq() += data.dImDdq.transpose() * s.beta + kkt_matrix_ref.Vq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += kkt_matrix_ref.Vv().transpose() * s.mu_stack(); 
  kkt_residual_ref.ldv  += data.dImDddv.transpose() * s.beta;
  kkt_residual_ref.lf() += dImDdf.transpose() * s.beta;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  const double l1norm_ref = data.ImD.lpNorm<1>() + kkt_residual_ref.V().lpNorm<1>();
  const double squarednorm_ref = data.ImD.squaredNorm() + kkt_residual_ref.V().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm_ref, l1norm);
  EXPECT_DOUBLE_EQ(squarednorm_ref, squarednorm);
}  


void ImpulseDynamicsBackwardEulerTest::testCondensing(Robot& robot, const ImpulseStatus& impulse_status) {
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
  kkt_matrix.Qdvdvff().diagonal().setRandom();
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseDynamicsBackwardEuler id(robot);
  ImpulseDynamicsBackwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  data.dImDdq.setRandom();
  data.dImDddv.setRandom();
  data.Minv.setRandom();
  data.Minv.template triangularView<Eigen::StrictlyLower>() 
      = data.Minv.transpose().template triangularView<Eigen::StrictlyLower>();
  data.Qdvq.setRandom();
  data.ImD.setRandom();
  data.Minv_ImD.setRandom();
  data.ldv.setRandom();
  data.Qdvf().setRandom();
  ImpulseDynamicsBackwardEulerData data_ref = data;
  ImpulseDynamicsBackwardEuler::condensing(robot, data, kkt_matrix, kkt_residual);
  const Eigen::MatrixXd Minv_dImDdq = data_ref.Minv * data_ref.dImDdq;
  const Eigen::MatrixXd Minv_dImDdf = - data_ref.Minv * kkt_matrix_ref.Vv().transpose();
  kkt_matrix_ref.Fvq() = Minv_dImDdq;
  kkt_matrix_ref.Fvf() = Minv_dImDdf; 
  data_ref.Minv_ImD = data_ref.Minv * data_ref.ImD; 
  data_ref.Qdvq   = - kkt_matrix_ref.Qdvdv() * Minv_dImDdq;
  data_ref.Qdvf() = - kkt_matrix_ref.Qdvdv() * Minv_dImDdf;
  data_ref.ldv    = - kkt_matrix_ref.Qdvdv() * data_ref.Minv_ImD;
  kkt_matrix_ref.Qqq()  -= Minv_dImDdq.transpose() * data_ref.Qdvq;
  kkt_matrix_ref.Qfq().transpose()  -= Minv_dImDdq.transpose() * data_ref.Qdvf();
  kkt_matrix_ref.Qff()  -= Minv_dImDdf.transpose() * data_ref.Qdvf();
  kkt_residual_ref.lq() -= Minv_dImDdq.transpose() * data_ref.ldv;
  kkt_residual_ref.lf() -= Minv_dImDdf.transpose() * data_ref.ldv;
  kkt_matrix_ref.Fvv() = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_residual_ref.Fv() -= data_ref.Minv_ImD;
  EXPECT_TRUE(data_ref.Minv_ImD.isApprox(data.Minv_ImD));
  EXPECT_TRUE(data_ref.Qdvq.isApprox(data.Qdvq));
  EXPECT_TRUE(data_ref.Qdvf().isApprox(data.Qdvf()));
  EXPECT_TRUE(data_ref.ldv.isApprox(data.ldv));
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
}


void ImpulseDynamicsBackwardEulerTest::testIntegration(Robot& robot, const ImpulseStatus& impulse_status) {
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
  kkt_matrix.Qdvdvff().diagonal().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseSplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  robot.updateKinematics(s.q, s.v);
  ImpulseDynamicsBackwardEuler id(robot), id_ref(robot);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix, kkt_residual);
  ImpulseDynamicsBackwardEulerData data_ref(robot);
  data_ref.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  id_ref.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsBackwardEuler::linearizeInverseImpulseDynamics(robot, impulse_status, s, data_ref);
  ImpulseDynamicsBackwardEuler::linearizeImpulseVelocityConstraint(robot, impulse_status, kkt_matrix_ref, kkt_residual_ref);
  robot.computeMinv(data_ref.dImDddv, data_ref.Minv);
  ImpulseDynamicsBackwardEuler::condensing(robot, data_ref, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitDirection d_ref = d;
  id.computeCondensedPrimalDirection(robot, kkt_matrix, d);
  d_ref.ddv()  = - data_ref.Minv_ImD;
  d_ref.ddv() -= kkt_matrix_ref.Fvq() * d.dq();
  d_ref.ddv() -= kkt_matrix_ref.Fvf() * d.df();
  EXPECT_TRUE(d.isApprox(d_ref));
  id.computeCondensedDualDirection(robot, d);
  data_ref.ldv += data_ref.Qdvq * d_ref.dq();
  data_ref.ldv += data_ref.Qdvf() * d_ref.df();
  data_ref.ldv += d_ref.dgmm();
  d_ref.dbeta() = - data_ref.Minv * data_ref.ldv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ImpulseDynamicsBackwardEulerTest::testComputeResidual(Robot& robot, const ImpulseStatus& impulse_status) {
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ImpulseDynamicsBackwardEuler id(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  robot.updateKinematics(s.q, s.v);
  id.computeImpulseDynamicsResidual(robot, impulse_status, s, kkt_residual);
  const double l1norm = id.l1NormImpulseDynamicsResidual(kkt_residual);
  const double squarednorm = id.squaredNormImpulseDynamicsResidual(kkt_residual);
  ImpulseDynamicsBackwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data.ImD);
  robot.updateKinematics(s.q, s.v);
  robot.computeImpulseVelocityResidual(impulse_status, kkt_residual_ref.V());
  double l1norm_ref = data.ImD.lpNorm<1>() + kkt_residual_ref.V().lpNorm<1>();
  double squarednorm_ref = data.ImD.squaredNorm() + kkt_residual_ref.V().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ImpulseDynamicsBackwardEulerTest, fixedBase) {
  std::vector<int> impulse_frames = {18};
  ImpulseStatus impulse_status(impulse_frames.size());
  for (int i=0; i<impulse_frames.size(); ++i) {
    impulse_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(fixed_base_urdf, impulse_frames);
  impulse_status.setImpulseStatus({false});
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
}


TEST_F(ImpulseDynamicsBackwardEulerTest, floatingBase) {
  std::vector<int> impulse_frames = {14, 24, 34, 44};
  ImpulseStatus impulse_status(impulse_frames.size());
  for (int i=0; i<impulse_frames.size(); ++i) {
    impulse_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(floating_base_urdf, impulse_frames);
  impulse_status.setImpulseStatus({false, false, false, false});
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
  std::random_device rnd;
  std::vector<bool> is_impulse_active;
  for (const auto frame : impulse_frames) {
    is_impulse_active.push_back(rnd()%2==0);
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  impulse_status.setImpulseStatus(is_impulse_active);
  testLinearizeInverseImpulseDynamics(robot, impulse_status);
  testLinearizeImpulseVelocityConstraints(robot, impulse_status);
  testLinearizeImpulseDynamics(robot, impulse_status);
  testCondensing(robot, impulse_status);
  testIntegration(robot, impulse_status);
  testComputeResidual(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}