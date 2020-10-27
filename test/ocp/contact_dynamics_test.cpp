#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

class ContactDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condenseContactDynamicsFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix_ref.Qxx() = kkt_matrix.Qxx();
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix_ref.Qaaff() = kkt_matrix.Qaaff();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  Eigen::MatrixXd Minv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.computeMinv(dID_da, Minv);
  Eigen::MatrixXd Minv_dID_dqv = Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv());
  Minv_dID_dqv.leftCols(robot.dimv()) = Minv * dID_dq;
  Minv_dID_dqv.rightCols(robot.dimv()) = Minv * dID_dv;
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  const Eigen::MatrixXd Qafqv_condensed = - kkt_matrix_ref.Qaaff() * Minv_dID_dqv;
  const Eigen::MatrixXd Qafu_condensed = - kkt_matrix_ref.Qaaff() * Minv * IO_mat;
  Eigen::VectorXd laf = kkt_residual.la;
  const Eigen::VectorXd laf_condensed = laf - kkt_matrix_ref.Qaaff() * Minv * ID;
  kkt_matrix_ref.Qxx() -= Minv_dID_dqv.transpose() * Qafqv_condensed;
  kkt_matrix_ref.Qxu() -= Minv_dID_dqv.transpose() * Qafu_condensed;
  kkt_matrix_ref.Quu() += IO_mat.transpose() * Minv * Qafu_condensed;
  Eigen::VectorXd Minv_ID = Eigen::VectorXd::Zero(robot.dimv());
  kkt_residual_ref.lx() -= Minv_dID_dqv.transpose() * laf_condensed;
  kkt_residual_ref.lu() += IO_mat.transpose() * Minv * laf_condensed;
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix_ref.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxu().isApprox(kkt_matrix_ref.Qxu()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix_ref.Quu()));
  EXPECT_TRUE(kkt_matrix.Qux().isZero());
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  Eigen::MatrixXd IOO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IOO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau_ * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxx = IOO_mat * Minv_dID_dqv;
  const Eigen::MatrixXd Fxu = IOO_mat * Minv * IO_mat;
  const Eigen::VectorXd Fx = IOO_mat * Minv * ID;
  kkt_matrix_ref.Fvq() = - Fxx.bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - Fxx.bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = Fxu.bottomRows(robot.dimv());
  kkt_residual_ref.Fv() = - Fx.tail(robot.dimv());
  EXPECT_TRUE(kkt_matrix.Fvq().isApprox(kkt_matrix_ref.Fvq()));
  EXPECT_TRUE(kkt_matrix.Fvv().isApprox(kkt_matrix_ref.Fvv()));
  EXPECT_TRUE(kkt_matrix.Fvu().isApprox(kkt_matrix_ref.Fvu()));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  cd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d.dgmm(), d);
  Eigen::VectorXd da_ref = - Minv * ID;
  da_ref -= Minv_dID_dqv * d.dx();
  da_ref += Minv * IO_mat * d.du();
  EXPECT_TRUE(da_ref.isApprox(d.da()));
  Eigen::VectorXd laf_tmp = laf_condensed;
  laf_tmp.head(robot.dimv()) += dtau_ * d.dgmm();
  laf_tmp += Qafu_condensed * d.du();
  laf_tmp += Qafqv_condensed * d.dx();
  const Eigen::VectorXd dbeta_ref = - Minv * laf_tmp / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {false};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ContactDynamics cd(robot);
  cd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(contact_status, s.f);
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  const double l1norm_ref = dtau_ * ID.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condenseContactDynamicsFixedBaseWithContact) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix_ref.Qxx() = kkt_matrix.Qxx();
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix_ref.Qaaff() = kkt_matrix.Qaaff();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  const int dimf = contact_status.dimf();
  const int dimv = robot.dimv();
  const int dimvf = dimv + dimf;
  Eigen::MatrixXd MJtJinv = Eigen::MatrixXd::Zero(dimvf, dimvf);
  robot.computeMJtJinv(contact_status, dID_da, dC_da, MJtJinv);
  Eigen::MatrixXd dIDC_dqv = Eigen::MatrixXd::Zero(dimvf, 2*dimv);
  dIDC_dqv.topLeftCorner(dimv, dimv) = dID_dq;
  dIDC_dqv.topRightCorner(dimv, dimv) = dID_dv;
  dIDC_dqv.bottomLeftCorner(dimf, dimv) = dC_dq;
  dIDC_dqv.bottomRightCorner(dimf, dimv) = dC_dv;
  const Eigen::MatrixXd MJtJinv_dIDC_dqv = MJtJinv * dIDC_dqv;
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimvf, robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  const Eigen::MatrixXd Qafqv_condensed = - kkt_matrix_ref.Qaaff() * MJtJinv_dIDC_dqv;
  const Eigen::MatrixXd Qafu_condensed = - kkt_matrix_ref.Qaaff() * MJtJinv * IO_mat;
  Eigen::VectorXd laf = Eigen::VectorXd::Zero(dimvf);
  laf.head(dimv) = kkt_residual.la;
  laf.tail(dimf) = - kkt_residual.lf();
  Eigen::VectorXd IDC = Eigen::VectorXd::Zero(dimvf);
  IDC.head(dimv) = ID;
  IDC.tail(dimf) = C;
  const Eigen::VectorXd laf_condensed = laf - kkt_matrix_ref.Qaaff() * MJtJinv * IDC;
  kkt_matrix_ref.Qxx() -= MJtJinv_dIDC_dqv.transpose() * Qafqv_condensed;
  kkt_matrix_ref.Qxu_full() -= MJtJinv_dIDC_dqv.transpose() * Qafu_condensed;
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * MJtJinv * Qafu_condensed;
  kkt_residual_ref.lx() -= MJtJinv_dIDC_dqv.transpose() * laf_condensed;
  kkt_residual_ref.lu() += IO_mat.transpose() * MJtJinv * laf_condensed;
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix_ref.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxu().isApprox(kkt_matrix_ref.Qxu()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix_ref.Quu()));
  EXPECT_TRUE(kkt_matrix.Qux().isZero());
  Eigen::MatrixXd IOO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IOO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau_ * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxx = IOO_mat * MJtJinv_dIDC_dqv;
  const Eigen::MatrixXd Fxu = IOO_mat * MJtJinv * IO_mat;
  const Eigen::VectorXd Fx = IOO_mat * MJtJinv * IDC;
  kkt_matrix_ref.Fvq() = - Fxx.bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - Fxx.bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = Fxu.bottomRows(robot.dimv());
  kkt_residual_ref.Fv() = - Fx.tail(robot.dimv());
  EXPECT_TRUE(kkt_matrix.Fvq().isApprox(kkt_matrix_ref.Fvq()));
  EXPECT_TRUE(kkt_matrix.Fvv().isApprox(kkt_matrix_ref.Fvv()));
  EXPECT_TRUE(kkt_matrix.Fvu().isApprox(kkt_matrix_ref.Fvu()));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  cd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d.dgmm(), d);
  Eigen::VectorXd damf_ref = - MJtJinv * IDC;
  damf_ref -= MJtJinv_dIDC_dqv * d.dx();
  damf_ref += MJtJinv * IO_mat * d.du();
  EXPECT_TRUE(damf_ref.head(dimv).isApprox(d.da()));
  EXPECT_TRUE((-1*damf_ref.tail(dimf)).isApprox(d.df()));
  Eigen::VectorXd laf_tmp = laf_condensed;
  laf_tmp.head(robot.dimv()) += dtau_ * d.dgmm();
  laf_tmp += Qafu_condensed * d.du();
  laf_tmp += Qafqv_condensed * d.dx();
  const Eigen::VectorXd dbetamu_ref = - MJtJinv * laf_tmp / dtau_;
  EXPECT_TRUE(dbetamu_ref.head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(dbetamu_ref.tail(dimf).isApprox(d.dmu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(contact_status, s.f);
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(false);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(robot.dim_passive()) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isZero());
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(kkt_residual_ref.lu_passive));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * ID.lpNorm<1>() + dtau_ * s.u_passive.lpNorm<1>();
  const double squarednorm_ref = dtau_ * dtau_ * ID.squaredNorm() + dtau_ * dtau_ * s.u_passive.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condenseContactDynamicsFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(false);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lu() = kkt_residual.lu();
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix_ref.Qxx() = kkt_matrix.Qxx();
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix_ref.Qaaff() = kkt_matrix.Qaaff();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(6) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  Eigen::MatrixXd Minv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.computeMinv(dID_da, Minv);
  Eigen::MatrixXd Minv_dID_dqv = Eigen::MatrixXd::Zero(robot.dimv(), 2*robot.dimv());
  Minv_dID_dqv.leftCols(robot.dimv()) = Minv * dID_dq;
  Minv_dID_dqv.rightCols(robot.dimv()) = Minv * dID_dv;
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(robot.dimv()+contact_status.dimf(), robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  const Eigen::MatrixXd Qafqv_condensed = - kkt_matrix_ref.Qaaff() * Minv_dID_dqv;
  const Eigen::MatrixXd Qafu_condensed = - kkt_matrix_ref.Qaaff() * Minv * IO_mat;
  Eigen::VectorXd laf = kkt_residual.la;
  const Eigen::VectorXd laf_condensed = laf - kkt_matrix_ref.Qaaff() * Minv * ID;
  kkt_matrix_ref.Qxx() -= Minv_dID_dqv.transpose() * Qafqv_condensed;
  kkt_matrix_ref.Qxu_full() -= Minv_dID_dqv.transpose() * Qafu_condensed;
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * Minv * Qafu_condensed;
  Eigen::VectorXd Minv_ID = Eigen::VectorXd::Zero(robot.dimv());
  kkt_residual_ref.lx() -= Minv_dID_dqv.transpose() * laf_condensed;
  kkt_residual_ref.lu_passive += (IO_mat.transpose() * Minv * laf_condensed).head(robot.dim_passive());
  kkt_residual_ref.lu() += (IO_mat.transpose() * Minv * laf_condensed).tail(robot.dimu());
  kkt_residual_ref.lu() -= kkt_matrix_ref.Quu_full().bottomLeftCorner(robot.dimu(), robot.dim_passive()) * s.u_passive;
  kkt_residual_ref.lx() -= kkt_matrix_ref.Qxu_full().leftCols(robot.dim_passive()) * s.u_passive;
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix_ref.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxu().isApprox(kkt_matrix_ref.Qxu()));
  EXPECT_TRUE(kkt_matrix.Qxu_full().isApprox(kkt_matrix_ref.Qxu_full()));
  EXPECT_TRUE(kkt_matrix.Quu_full().isApprox(kkt_matrix_ref.Quu_full()));
  EXPECT_TRUE(kkt_matrix.Qux_full().isZero());
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(kkt_residual_ref.lu_passive));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  Eigen::MatrixXd IOO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IOO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau_ * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxx = IOO_mat * Minv_dID_dqv;
  const Eigen::MatrixXd Fxu_full = IOO_mat * Minv * IO_mat;
  const Eigen::VectorXd Fx = IOO_mat * Minv * ID;
  kkt_matrix_ref.Fvq() = - Fxx.bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - Fxx.bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = Fxu_full.bottomRightCorner(robot.dimv(), robot.dimu());
  kkt_residual_ref.Fv() -= Fx.tail(robot.dimv());
  kkt_residual_ref.Fv() -= Fxu_full.bottomLeftCorner(robot.dimv(), 6) * s.u_passive;
  EXPECT_TRUE(kkt_matrix.Fvq().isApprox(kkt_matrix_ref.Fvq()));
  EXPECT_TRUE(kkt_matrix.Fvv().isApprox(kkt_matrix_ref.Fvv()));
  EXPECT_TRUE(kkt_matrix.Fvu().isApprox(kkt_matrix_ref.Fvu()));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  cd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d.dgmm(), d);
  const Eigen::VectorXd du_passive_ref = - s.u_passive;
  Eigen::VectorXd dlmdgmm = Eigen::VectorXd::Zero(2*robot.dimv());
  dlmdgmm.head(robot.dimv()) = d.dlmd();
  dlmdgmm.tail(robot.dimv()) = d.dgmm();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(robot.dimv());
  du_full.head(robot.dim_passive()) = du_passive_ref;
  du_full.tail(robot.dimu()) = d.du();
  const Eigen::VectorXd m_dnu_passive_dtau
      = kkt_residual_ref.lu_passive 
        + kkt_matrix_ref.Qxu_full().transpose().topRows(robot.dim_passive()) * d.dx()
        + kkt_matrix_ref.Quu_full().topRows(robot.dim_passive()) * du_full 
        + IO_mat.leftCols(robot.dim_passive()).transpose() * Minv * IOO_mat.transpose() * dlmdgmm;
  const Eigen::VectorXd dnu_passive_ref = - m_dnu_passive_dtau / dtau_;
  EXPECT_TRUE(du_passive_ref.isApprox(d.du_passive));
  EXPECT_TRUE(dnu_passive_ref.isApprox(d.dnu_passive));
  Eigen::VectorXd da_ref = - Minv * ID;
  da_ref -= Minv_dID_dqv * d.dx();
  da_ref += Minv * IO_mat * du_full;
  EXPECT_TRUE(da_ref.isApprox(d.da()));
  Eigen::VectorXd laf_tmp = laf_condensed;
  laf_tmp.head(robot.dimv()) += dtau_ * d.dgmm();
  laf_tmp += Qafu_condensed * du_full;
  laf_tmp += Qafqv_condensed * d.dx();
  const Eigen::VectorXd dbeta_ref = - Minv * laf_tmp / dtau_;
  EXPECT_TRUE(dbeta_ref.isApprox(d.dbeta()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(false);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(contact_status, s.f);
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(6) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, linearizeContactDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(robot.dim_passive()) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.lu_passive.isApprox(kkt_residual_ref.lu_passive));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + s.u_passive.lpNorm<1>() + C.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + s.u_passive.squaredNorm() + C.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, condenseContactDynamicsFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual.lq() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lv() = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.la = Eigen::VectorXd::Random(robot.dimv());
  kkt_residual.lf() = Eigen::VectorXd::Random(contact_status.dimf());
  kkt_residual.lu() = Eigen::VectorXd::Random(robot.dimu());
  kkt_residual_ref.lq() = kkt_residual.lq();
  kkt_residual_ref.lv() = kkt_residual.lv();
  kkt_residual_ref.la = kkt_residual.la;
  kkt_residual_ref.lf() = kkt_residual.lf();
  kkt_residual_ref.lu() = kkt_residual.lu();
  kkt_matrix.Qxx() = Eigen::MatrixXd::Random(2*robot.dimv(), 2*robot.dimv());
  kkt_matrix_ref.Qxx() = kkt_matrix.Qxx();
  kkt_matrix.Qaaff().diagonal() = Eigen::VectorXd::Random(robot.dimv()+contact_status.dimf());
  kkt_matrix_ref.Qaaff() = kkt_matrix.Qaaff();
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.condenseContactDynamics(robot, contact_status, dtau_, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(robot.dim_passive()) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual_ref.lq() += dtau_ * dID_dq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau_ * dID_dv.transpose() * s.beta; 
  kkt_residual_ref.la += dtau_ * dID_da.transpose() * s.beta;
  kkt_residual_ref.lf() += dtau_ * dID_df.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau_ * s.beta.tail(robot.dimu());
  kkt_residual_ref.lu_passive -= dtau_ * s.beta.head(robot.dim_passive());
  kkt_residual_ref.lu_passive += dtau_ * s.nu_passive;
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau_, dC_dq, dC_dv, dC_da);
  kkt_residual_ref.lq() += dtau_ * dC_dq.transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau_ * dC_dv.transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau_ * dC_da.transpose() * s.mu_stack();
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimvf = dimv + dimf;
  Eigen::MatrixXd MJtJinv = Eigen::MatrixXd::Zero(dimvf, dimvf);
  robot.computeMJtJinv(contact_status, dID_da, dC_da, MJtJinv);
  Eigen::MatrixXd dIDC_dqv = Eigen::MatrixXd::Zero(dimvf, 2*dimv);
  dIDC_dqv.topLeftCorner(dimv, dimv) = dID_dq;
  dIDC_dqv.topRightCorner(dimv, dimv) = dID_dv;
  dIDC_dqv.bottomLeftCorner(dimf, dimv) = dC_dq;
  dIDC_dqv.bottomRightCorner(dimf, dimv) = dC_dv;
  const Eigen::MatrixXd MJtJinv_dIDC_dqv = MJtJinv * dIDC_dqv;
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimvf, robot.dimv());
  IO_mat.topRows(robot.dimv()).setIdentity();
  const Eigen::MatrixXd Qafqv_condensed = - kkt_matrix_ref.Qaaff() * MJtJinv_dIDC_dqv;
  const Eigen::MatrixXd Qafu_condensed = - kkt_matrix_ref.Qaaff() * MJtJinv * IO_mat;
  Eigen::VectorXd laf = Eigen::VectorXd::Zero(dimvf);
  laf.head(dimv) = kkt_residual.la;
  laf.tail(dimf) = - kkt_residual.lf();
  Eigen::VectorXd IDC = Eigen::VectorXd::Zero(dimvf);
  IDC.head(dimv) = ID;
  IDC.tail(dimf) = C;
  const Eigen::VectorXd laf_condensed = laf - kkt_matrix_ref.Qaaff() * MJtJinv * IDC;
  kkt_matrix_ref.Qxx() -= MJtJinv_dIDC_dqv.transpose() * Qafqv_condensed;
  kkt_matrix_ref.Qxu_full() -= MJtJinv_dIDC_dqv.transpose() * Qafu_condensed;
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * MJtJinv * Qafu_condensed;
  kkt_residual_ref.lx() -= MJtJinv_dIDC_dqv.transpose() * laf_condensed;
  kkt_residual_ref.lu_passive += (IO_mat.transpose() * MJtJinv * laf_condensed).head(robot.dim_passive());
  kkt_residual_ref.lu() += (IO_mat.transpose() * MJtJinv * laf_condensed).tail(robot.dimu());
  kkt_residual_ref.lu() -= kkt_matrix_ref.Quu_full().bottomLeftCorner(robot.dimu(), robot.dim_passive()) * s.u_passive;
  kkt_residual_ref.lx() -= kkt_matrix_ref.Qxu_full().leftCols(robot.dim_passive()) * s.u_passive;
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix_ref.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxu().isApprox(kkt_matrix_ref.Qxu()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix_ref.Quu()));
  EXPECT_TRUE(kkt_matrix.Qux().isZero());
  Eigen::MatrixXd IOO_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dimv()+contact_status.dimf());
  IOO_mat.bottomLeftCorner(robot.dimv(), robot.dimv()) = dtau_ * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  const Eigen::MatrixXd Fxx = IOO_mat * MJtJinv_dIDC_dqv;
  const Eigen::MatrixXd Fxu_full = IOO_mat * MJtJinv * IO_mat;
  const Eigen::VectorXd Fx = IOO_mat * MJtJinv * IDC;
  kkt_matrix_ref.Fvq() = - Fxx.bottomLeftCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv()) - Fxx.bottomRightCorner(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fvu() = Fxu_full.bottomRightCorner(robot.dimv(), robot.dimu());
  kkt_residual_ref.Fv() -= Fx.tail(robot.dimv());
  kkt_residual_ref.Fv() -= Fxu_full.bottomLeftCorner(robot.dimv(), 6) * s.u_passive;
  EXPECT_TRUE(kkt_matrix.Fvq().isApprox(kkt_matrix_ref.Fvq()));
  EXPECT_TRUE(kkt_matrix.Fvv().isApprox(kkt_matrix_ref.Fvv()));
  EXPECT_TRUE(kkt_matrix.Fvu().isApprox(kkt_matrix_ref.Fvu()));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isApprox(kkt_residual_ref.Fv()));
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  cd.computeCondensedDirection(dtau_, kkt_matrix, kkt_residual, d.dgmm(), d);
  const Eigen::VectorXd du_passive_ref = - s.u_passive;
  Eigen::VectorXd dlmdgmm = Eigen::VectorXd::Zero(2*robot.dimv());
  dlmdgmm.head(robot.dimv()) = d.dlmd();
  dlmdgmm.tail(robot.dimv()) = d.dgmm();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(robot.dimv());
  du_full.head(robot.dim_passive()) = du_passive_ref;
  du_full.tail(robot.dimu()) = d.du();
  const Eigen::VectorXd m_dnu_passive_dtau
      = kkt_residual_ref.lu_passive 
        + kkt_matrix_ref.Qxu_full().transpose().topRows(robot.dim_passive()) * d.dx()
        + kkt_matrix_ref.Quu_full().topRows(robot.dim_passive()) * du_full 
        + IO_mat.leftCols(robot.dim_passive()).transpose() * MJtJinv * IOO_mat.transpose() * dlmdgmm;
  const Eigen::VectorXd dnu_passive_ref = - m_dnu_passive_dtau / dtau_;
  EXPECT_TRUE(du_passive_ref.isApprox(d.du_passive));
  EXPECT_TRUE(dnu_passive_ref.isApprox(d.dnu_passive));
  Eigen::VectorXd damf_ref = - MJtJinv * IDC;
  damf_ref -= MJtJinv_dIDC_dqv * d.dx();
  damf_ref += MJtJinv * IO_mat * du_full;
  EXPECT_TRUE(damf_ref.head(dimv).isApprox(d.da()));
  EXPECT_TRUE((-1*damf_ref.tail(dimf)).isApprox(d.df()));
  Eigen::VectorXd laf_tmp = laf_condensed;
  laf_tmp.head(robot.dimv()) += dtau_ * d.dgmm();
  laf_tmp += Qafu_condensed * du_full;
  laf_tmp += Qafqv_condensed * d.dx();
  const Eigen::VectorXd dbetamu_ref = - MJtJinv * laf_tmp / dtau_;
  EXPECT_TRUE(dbetamu_ref.head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(dbetamu_ref.tail(dimf).isApprox(d.dmu()));
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, computeContactDynamicsResidualFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  SplitSolution s = SplitSolution::Random(robot, contact_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  KKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  KKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.computeContactDynamicsResidual(robot, contact_status, dtau_, s, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau_);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau_);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.setContactForces(contact_status, s.f);
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.head(6) -= s.u_passive;
  ID.tail(robot.dimu()) -= s.u;
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.computeBaumgarteResidual(contact_status, dtau_, C);
  const double l1norm_ref = dtau_ * (ID.lpNorm<1>() + C.lpNorm<1>() + s.u_passive.lpNorm<1>());
  const double squarednorm_ref = dtau_ * dtau_ * (ID.squaredNorm() + C.squaredNorm() + s.u_passive.squaredNorm());
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}