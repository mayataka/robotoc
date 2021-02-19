#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/contact_dynamics_data.hpp"
#include "idocp/ocp/contact_dynamics.hpp"


namespace idocp {

class ContactDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    baumgarte_time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testLinearizeInverseDynamics(Robot& robot, const ContactStatus& contact_status) const;
  void testLinearizeContactConstraints(Robot& robot, const ContactStatus& contact_status) const;
  void testLinearizeContactDynamics(Robot& robot, const ContactStatus& contact_status) const;
  void testCondenseContactDynamics(Robot& robot, const ContactStatus& contact_status, const bool is_forward_euler) const;
  void testComputeCondensedPrimalDirection(Robot& robot, const ContactStatus& contact_status) const;
  void testComputeCondensedDualDirection(Robot& robot, const ContactStatus& contact_status) const;
  void testComputeResidual(Robot& robot, const ContactStatus& contact_status) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double baumgarte_time_step, dtau;
};


void ContactDynamicsTest::testLinearizeInverseDynamics(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID().noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


void ContactDynamicsTest::testLinearizeContactConstraints(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  data.setContactStatus(contact_status);
  data_ref.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, baumgarte_time_step, data);
  if (contact_status.hasActiveContacts()) {
    robot.computeBaumgarteResidual(contact_status, baumgarte_time_step, contact_status.contactPoints(), data_ref.C());
    robot.computeBaumgarteDerivatives(contact_status, baumgarte_time_step, data_ref.dCdq(), 
                                      data_ref.dCdv(), data_ref.dCda());
    EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
    EXPECT_TRUE(data.dCda().isApprox(data_ref.dCda()));
    EXPECT_TRUE(data.dIDCdqv().isApprox(data_ref.dIDCdqv()));
  }
  else {
    EXPECT_TRUE(data.IDC().isZero());
    EXPECT_TRUE(data.dCda().isZero());
    EXPECT_TRUE(data.dIDCdqv().isZero());
  }
}


void ContactDynamicsTest::testLinearizeContactDynamics(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom();
  kkt_residual.la.setRandom();
  kkt_residual.lf().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.lu_passive.setRandom();
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(robot.dimv());
  lu_full_ref.head(robot.dim_passive()).setZero();
  lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, baumgarte_time_step, data);
  Eigen::MatrixXd dIDdf = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dIDdf);
  kkt_residual_ref.lq() += dtau * data.dIDdq().transpose() * s.beta + dtau * data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += dtau * data.dIDdv().transpose() * s.beta + dtau * data.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.la += dtau * data.dIDda.transpose() * s.beta + dtau * data.dCda().transpose() * s.mu_stack();
  if (contact_status.hasActiveContacts()) {
    kkt_residual_ref.lf() += dtau * dIDdf.transpose() * s.beta;
  }
  lu_full_ref -= dtau * s.beta;
  if (robot.hasFloatingBase()) {
    lu_full_ref.head(robot.dim_passive()) += dtau * s.nu_passive;
  }
  EXPECT_TRUE(kkt_residual.lq().isApprox(kkt_residual_ref.lq()));
  EXPECT_TRUE(kkt_residual.lv().isApprox(kkt_residual_ref.lv()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  if (contact_status.hasActiveContacts()) {
    EXPECT_TRUE(kkt_residual.lf().isApprox(kkt_residual_ref.lf()));
  }
  if (robot.hasFloatingBase()) {
    EXPECT_TRUE(kkt_residual.lu_passive.isApprox(lu_full_ref.head(robot.dim_passive())));
    EXPECT_TRUE(kkt_residual.lu().isApprox(lu_full_ref.tail(robot.dimu())));
  }
  else {
    EXPECT_TRUE(kkt_residual.lu().isApprox(lu_full_ref));
  }
  const double l1norm_ref = dtau * data.IDC().lpNorm<1>();
  const double squarednorm_ref = dtau * dtau * data.IDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}  


void ContactDynamicsTest::testCondenseContactDynamics(Robot& robot, const ContactStatus& contact_status, const bool is_forward_euler) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lx().setRandom();
  kkt_residual.la.setRandom();
  kkt_residual.lf().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.setZero();
  }
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qxu_full().setRandom();
  kkt_matrix.Qux_full() = kkt_matrix.Qxu_full().transpose();
  kkt_matrix.Quu_full().setRandom();
  kkt_matrix.Quu_full().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Quu_full().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qaa().diagonal().setRandom();
  const Eigen::MatrixXd Qff_seed = Eigen::MatrixXd::Random(dimf, dimf);
  kkt_matrix.Qff() = Qff_seed * Qff_seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(dim_passive, dim_passive).setRandom();
  }
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix, kkt_residual, is_forward_euler);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, baumgarte_time_step, data);
  robot.computeMJtJinv(data.dIDda, data.dCda(), data.MJtJinv());
  data.MJtJinv_dIDCdqv() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC() = data.MJtJinv() * data.IDC();
  data.Qafqv() = - kkt_matrix_ref.Qaaff() * data.MJtJinv_dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  data.Qafu_full() = kkt_matrix_ref.Qaaff() * data.MJtJinv() * IO_mat;
  data.la() = kkt_residual_ref.la;
  data.lf() = - kkt_residual_ref.lf();
  data.laf() -= kkt_matrix_ref.Qaaff() * data.MJtJinv() * data.IDC();
  kkt_matrix_ref.Qxx() -= data.MJtJinv_dIDCdqv().transpose() * data.Qafqv();
  kkt_matrix_ref.Qxu_full() -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full();
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * data.MJtJinv() * data.Qafu_full();
  kkt_residual_ref.lx() -= data.MJtJinv_dIDCdqv().transpose() * data.laf();
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(dimv);
  if (robot.hasFloatingBase()) {
    lu_full_ref.head(dim_passive) = kkt_residual_ref.lu_passive;
    lu_full_ref.tail(dimu) = kkt_residual_ref.lu();
  }
  else {
    lu_full_ref = kkt_residual_ref.lu();
  }
  lu_full_ref += IO_mat.transpose() * data.MJtJinv() * data.laf();
  if (robot.hasFloatingBase()) {
    kkt_residual_ref.lu_passive = lu_full_ref.head(dim_passive);
    kkt_residual_ref.lu() = lu_full_ref.tail(dimu);
  }
  else {
    kkt_residual_ref.lu() = lu_full_ref;
  }
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  if (is_forward_euler) {
    kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv);
  }
  else {
    kkt_matrix_ref.Fvv() = - Eigen::MatrixXd::Identity(dimv, dimv);
  }
  kkt_matrix_ref.Fxx() -= OOIO_mat * data.MJtJinv_dIDCdqv();
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fxu() = Fxu_full.rightCols(dimu);
  kkt_residual_ref.Fx() -= (OOIO_mat * data.MJtJinv() * data.IDC());
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
}


void ContactDynamicsTest::testComputeCondensedPrimalDirection(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lx().setRandom();
  kkt_residual.la.setRandom();
  kkt_residual.lf().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.setZero();
  }
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qxu_full().setRandom();
  kkt_matrix.Qux_full() = kkt_matrix.Qxu_full().transpose();
  kkt_matrix.Quu_full().setRandom();
  kkt_matrix.Quu_full().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Quu_full().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qaa().diagonal().setRandom();
  const Eigen::MatrixXd Qff_seed = Eigen::MatrixXd::Random(dimf, dimf);
  kkt_matrix.Qff() = Qff_seed * Qff_seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(dim_passive, dim_passive).setRandom();
  }
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix, kkt_residual, true);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, baumgarte_time_step, data);
  robot.computeMJtJinv(data.dIDda, data.dCda(), data.MJtJinv());
  data.MJtJinv_dIDCdqv() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC() = data.MJtJinv() * data.IDC();
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  cd.computeCondensedPrimalDirection(robot, d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  du_full.tail(robot.dimu()) = d_ref.du();
  d_ref.daf() = - data.MJtJinv() * (data.dIDCdqv() * d.dx() - IO_mat * du_full + data.IDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ContactDynamicsTest::testComputeCondensedDualDirection(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.lx().setRandom();
  kkt_residual.la.setRandom();
  kkt_residual.lf().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.setZero();
  }
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qxu_full().setRandom();
  kkt_matrix.Qux_full() = kkt_matrix.Qxu_full().transpose();
  kkt_matrix.Quu_full().setRandom();
  kkt_matrix.Quu_full().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Quu_full().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qaa().diagonal().setRandom();
  const Eigen::MatrixXd Qff_seed = Eigen::MatrixXd::Random(dimf, dimf);
  kkt_matrix.Qff() = Qff_seed * Qff_seed.transpose();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(dim_passive, dim_passive).setRandom();
  }
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix, kkt_residual, true);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, baumgarte_time_step, data);
  robot.computeMJtJinv(data.dIDda, data.dCda(), data.MJtJinv());
  data.MJtJinv_dIDCdqv() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC() = data.MJtJinv() * data.IDC();
  data.Qafqv() = - kkt_matrix.Qaaff() * data.MJtJinv_dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  data.Qafu_full() = kkt_matrix.Qaaff() * data.MJtJinv() * IO_mat;
  data.la() = kkt_residual.la;
  data.lf() = - kkt_residual.lf();
  data.laf() -= kkt_matrix.Qaaff() * data.MJtJinv() * data.IDC();
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const Eigen::VectorXd dlmdgmm = Eigen::VectorXd::Random(2*dimv);
  const Eigen::VectorXd dgmm = dlmdgmm.tail(dimv);
  cd.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, dgmm, d);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  if (robot.hasFloatingBase()) {
    du_full.tail(dimu) = d_ref.du(); 
    d_ref.dnu_passive = - (kkt_residual.lu_passive + (kkt_matrix.Qxu_full().leftCols(dim_passive)).transpose() * d_ref.dx() 
                            + kkt_matrix.Quu_full().topRows(dim_passive) * du_full
                            + (IO_mat.transpose() * data.MJtJinv() * OOIO_mat.transpose() * dlmdgmm).head(dim_passive)) / dtau;
  }
  else {
    du_full = d_ref.du(); 
  }
  d_ref.dbetamu() = - data.MJtJinv() * (data.Qafqv() * d_ref.dx() 
                                        + data.Qafu_full() * du_full 
                                        + OOIO_mat.transpose() * dlmdgmm
                                        + data.laf()) / dtau;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ContactDynamicsTest::testComputeResidual(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.computeContactDynamicsResidual(robot, contact_status, s);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  if (contact_status.hasActiveContacts()) {
    robot.computeBaumgarteResidual(contact_status, baumgarte_time_step, contact_status.contactPoints(), data.C());
  }
  double l1norm_ref = dtau * data.IDC().lpNorm<1>();
  double squarednorm_ref = dtau * dtau * data.IDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondenseContactDynamics(robot, contact_status, true);
  testCondenseContactDynamics(robot, contact_status, false);
  testComputeCondensedPrimalDirection(robot, contact_status);
  testComputeCondensedDualDirection(robot, contact_status);
  testComputeResidual(robot, contact_status);
  contact_status.activateContact(0);
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testCondenseContactDynamics(robot, contact_status, true);
  testCondenseContactDynamics(robot, contact_status, false);
  testComputeCondensedPrimalDirection(robot, contact_status);
  testComputeCondensedDualDirection(robot, contact_status);
  testComputeResidual(robot, contact_status);
}


TEST_F(ContactDynamicsTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondenseContactDynamics(robot, contact_status, true);
  testCondenseContactDynamics(robot, contact_status, false);
  testComputeCondensedPrimalDirection(robot, contact_status);
  testComputeCondensedDualDirection(robot, contact_status);
  testComputeResidual(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondenseContactDynamics(robot, contact_status, true);
  testCondenseContactDynamics(robot, contact_status, false);
  testComputeCondensedPrimalDirection(robot, contact_status);
  testComputeCondensedDualDirection(robot, contact_status);
  testComputeResidual(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}