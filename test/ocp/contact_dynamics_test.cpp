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
  }

  virtual void TearDown() {
  }

  static void testLinearizeInverseDynamics(Robot& robot, const ContactStatus& contact_status);
  static void testLinearizeContactConstraints(Robot& robot, const ContactStatus& contact_status);
  static void testLinearizeContactDynamics(Robot& robot, const ContactStatus& contact_status);
  static void testCondensing(Robot& robot, const ContactStatus& contact_status);
  static void testExpansionPrimal(Robot& robot, const ContactStatus& contact_status);
  static void testExpansionDual(Robot& robot, const ContactStatus& contact_status);
  static void testIntegration(Robot& robot, const ContactStatus& contact_status);
  static void testComputeResidual(Robot& robot, const ContactStatus& contact_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void ContactDynamicsTest::testLinearizeInverseDynamics(Robot& robot, const ContactStatus& contact_status) {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID().noalias() -= s.u;
  if (robot.hasFloatingBase()) {
    data_ref.ID_passive().noalias() -= s.u_passive;
  }
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  EXPECT_TRUE(data_ref.IDC().isApprox(data.IDC()));
  EXPECT_TRUE(data_ref.dIDda.isApprox(data.dIDda));
  EXPECT_TRUE(data_ref.dIDdq().isApprox(data.dIDdq()));
  EXPECT_TRUE(data_ref.dIDdv().isApprox(data.dIDdv()));
}


void ContactDynamicsTest::testLinearizeContactConstraints(Robot& robot, const ContactStatus& contact_status) {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot), data_ref(robot);
  data.setContactStatus(contact_status);
  data_ref.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
  if (contact_status.hasActiveContacts()) {
    robot.computeBaumgarteResidual(contact_status, dtau, contact_status.contactPoints(), data_ref.C());
    robot.computeBaumgarteDerivatives(contact_status, dtau, data_ref.dCdq(), 
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


void ContactDynamicsTest::testLinearizeContactDynamics(Robot& robot, const ContactStatus& contact_status) {
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
  if (robot.hasFloatingBase()) {
    lu_full_ref.head(robot.dim_passive()) = kkt_residual_ref.lu_passive;
    lu_full_ref.tail(robot.dimu()) = kkt_residual_ref.lu();
  }
  else {
    lu_full_ref = kkt_residual_ref.lu();
  }
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data);
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
    const double l1norm_ref = dtau * (data.IDC().lpNorm<1>() + s.u_passive.lpNorm<1>());
    const double squarednorm_ref = dtau * dtau * (data.IDC().squaredNorm() + s.u_passive.squaredNorm());
    EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
    EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
  }
  else {
    EXPECT_TRUE(kkt_residual.lu().isApprox(lu_full_ref));
    const double l1norm_ref = dtau * (data.IDC().lpNorm<1>());
    const double squarednorm_ref = dtau * dtau * (data.IDC().squaredNorm());
    EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
    EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
  }
}  


void ContactDynamicsTest::testCondensing(Robot& robot, const ContactStatus& contact_status) {
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
    kkt_residual.lu_passive.setRandom();
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
  kkt_matrix.Qaaff().diagonal().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(dim_passive, dim_passive).setRandom();
  }
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  ContactDynamics cd(robot);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  if (robot.hasFloatingBase()) {
    data.u_passive.setRandom();
  }
  data.dCda().setRandom();
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.MJtJinv_dIDCdqv().setRandom();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.IDC().setRandom();
  data.MJtJinv_IDC().setRandom();
  data.laf().setRandom();
  ContactDynamicsData data_ref = data;
  data_ref.MJtJinv_dIDCdqv().setZero();
  data_ref.Qafqv().setZero();
  data_ref.Qafu_full().setZero();
  data_ref.MJtJinv_IDC().setZero();
  data_ref.laf().setZero();
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ContactDynamics::condensing(robot, dtau, data, kkt_matrix, kkt_residual);
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC() = data_ref.MJtJinv() * data_ref.IDC();
  data_ref.Qafqv() = - kkt_matrix_ref.Qaaff() * data_ref.MJtJinv_dIDCdqv();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  data_ref.Qafu_full() = kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * IO_mat;
  data_ref.laf().head(dimv) = kkt_residual_ref.la;
  data_ref.laf().tail(dimf) = - kkt_residual_ref.lf();
  data_ref.laf() -= kkt_matrix_ref.Qaaff() * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxu_full() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu_full();
  kkt_matrix_ref.Quu_full() += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu_full();
  kkt_residual_ref.lx() -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  Eigen::VectorXd lu_full_ref = Eigen::VectorXd::Zero(dimv);
  if (robot.hasFloatingBase()) {
    lu_full_ref.head(dim_passive) = kkt_residual_ref.lu_passive;
    lu_full_ref.tail(dimu) = kkt_residual_ref.lu();
  }
  else {
    lu_full_ref = kkt_residual_ref.lu();
  }
  lu_full_ref += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  if (robot.hasFloatingBase()) {
    kkt_residual_ref.lu_passive = lu_full_ref.head(dim_passive);
    kkt_residual_ref.lu() = lu_full_ref.tail(dimu);
    kkt_residual_ref.lu() -= kkt_matrix_ref.Quu_passive_bottomLeft() * data_ref.u_passive;
    kkt_residual_ref.lx() -= kkt_matrix_ref.Qxu_passive() * data_ref.u_passive;
  }
  else {
    kkt_residual_ref.lu() = lu_full_ref;
  }
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix_ref.Fvv().setIdentity();
  kkt_matrix_ref.Fxx() -= OOIO_mat * data_ref.MJtJinv_dIDCdqv();
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data_ref.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fxu() = Fxu_full.rightCols(dimu);
  kkt_residual_ref.Fx() -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());
  if (robot.hasFloatingBase()) {
    kkt_residual_ref.Fx() -= Fxu_full.leftCols(dim_passive) * data_ref.u_passive;
  }
  EXPECT_TRUE(data_ref.MJtJinv().isApprox(data.MJtJinv()));
  EXPECT_TRUE(data_ref.MJtJinv_dIDCdqv().isApprox(data.MJtJinv_dIDCdqv()));
  EXPECT_TRUE(data_ref.MJtJinv_IDC().isApprox(data.MJtJinv_IDC()));
  EXPECT_TRUE(data_ref.Qafqv().isApprox(data.Qafqv()));
  EXPECT_TRUE(data_ref.Qafu_full().isApprox(data.Qafu_full()));
  EXPECT_TRUE(data_ref.laf().isApprox(data.laf()));
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
}


void ContactDynamicsTest::testExpansionPrimal(Robot& robot, const ContactStatus& contact_status) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  if (robot.hasFloatingBase()) {
    data.u_passive.setRandom();
  }
  data.dIDCdqv().setRandom();
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.IDC().setRandom();
  data.MJtJinv_dIDCdqv() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC() = data.MJtJinv() * data.IDC();
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  ContactDynamics::expansionPrimal(robot, data, d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  if (robot.hasFloatingBase()) {
    d_ref.du_passive = - data.u_passive;
    du_full.head(dim_passive) = d_ref.du_passive; 
    du_full.tail(dimu) = d_ref.du(); 
  }
  else {
    du_full = d_ref.du(); 
  }
  d_ref.daf() = - data.MJtJinv() * (data.dIDCdqv() * d.dx() - IO_mat * du_full + data.IDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ContactDynamicsTest::testExpansionDual(Robot& robot, const ContactStatus& contact_status) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.setRandom();
  }
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_matrix.Qxu_full().setRandom();
  kkt_matrix.Quu_full().setRandom();
  kkt_matrix.Quu_full().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Quu_full().transpose().template triangularView<Eigen::StrictlyLower>();
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  data.dIDda.setRandom();
  if (robot.hasFloatingBase()) {
    data.u_passive.setRandom();
  }
  data.MJtJinv().setRandom();
  data.MJtJinv().template triangularView<Eigen::StrictlyLower>() 
      = data.MJtJinv().transpose().template triangularView<Eigen::StrictlyLower>();
  data.Qafqv().setRandom();
  data.Qafu_full().setRandom();
  data.laf().setRandom();
  const ContactDynamicsData data_ref = data;
  const Eigen::VectorXd dlmdgmm_next = Eigen::VectorXd::Random(2*dimv);
  const Eigen::VectorXd dgmm_next = dlmdgmm_next.tail(dimv);
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ContactDynamics::expansionDual(robot, dtau, data, kkt_matrix, kkt_residual, dgmm_next, d);
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  if (robot.hasFloatingBase()) {
    du_full.head(dim_passive) = d_ref.du_passive; 
    du_full.tail(dimu) = d_ref.du(); 
    d_ref.dnu_passive = - (kkt_residual.lu_passive + (kkt_matrix.Qxu_full().leftCols(dim_passive)).transpose() * d_ref.dx() 
                            + kkt_matrix.Quu_full().topRows(dim_passive) * du_full
                            + (IO_mat.transpose() * data_ref.MJtJinv() * OOIO_mat.transpose() *dlmdgmm_next).head(dim_passive)) / dtau;
  }
  else {
    du_full = d_ref.du(); 
  }
  d_ref.dbetamu() = - data_ref.MJtJinv() * (data_ref.Qafqv() * d.dx() 
                                            + data_ref.Qafu_full() * du_full 
                                            + OOIO_mat.transpose() * dlmdgmm_next 
                                            + data_ref.laf()) / dtau;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ContactDynamicsTest::testIntegration(Robot& robot, const ContactStatus& contact_status) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.lx().setRandom();
  kkt_residual.la.setRandom();
  kkt_residual.lu().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.setRandom();
  }
  kkt_residual.Fx().setRandom();
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setContactStatus(contact_status);
  kkt_matrix.Qxx().setRandom();
  kkt_matrix.Qxx().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Qxx().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qxu_full().setRandom();
  kkt_matrix.Qux_full() = kkt_matrix.Qxu_full().transpose();
  kkt_matrix.Quu_full().setRandom();
  kkt_matrix.Quu_full().template triangularView<Eigen::StrictlyLower>()
      = kkt_matrix.Quu_full().transpose().template triangularView<Eigen::StrictlyLower>();
  kkt_matrix.Qaaff().diagonal().setRandom();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(dim_passive, dim_passive).setRandom();
  }
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics cd(robot), cd_ref(robot);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual);
  cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix, kkt_residual);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  if (robot.hasFloatingBase()) {
    data_ref.u_passive = s.u_passive;
  }
  cd_ref.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual_ref);
  ContactDynamics::linearizeInverseDynamics(robot, contact_status, s, data_ref);
  ContactDynamics::linearizeContactConstraint(robot, contact_status, dtau, data_ref);
  robot.computeMJtJinv(data_ref.dIDda, data_ref.dCda(), data_ref.MJtJinv());
  ContactDynamics::condensing(robot, dtau, data_ref, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitDirection d_ref = d;
  cd.computeCondensedPrimalDirection(robot, d);
  ContactDynamics::expansionPrimal(robot, data_ref, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  const Eigen::VectorXd dgmm_next = Eigen::VectorXd::Random(dimv);
  cd.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, dgmm_next, d);
  ContactDynamics::expansionDual(robot, dtau, data_ref, kkt_matrix_ref, kkt_residual_ref, dgmm_next, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ContactDynamicsTest::testComputeResidual(Robot& robot, const ContactStatus& contact_status) {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  cd.computeContactDynamicsResidual(robot, contact_status, dtau, s);
  const double l1norm = cd.l1NormContactDynamicsResidual(dtau);
  const double squarednorm = cd.squaredNormContactDynamicsResidual(dtau);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  if (robot.hasFloatingBase()) {
    data.ID_passive().noalias() -= s.u_passive;
  }
  robot.updateKinematics(s.q, s.v, s.a);
  if (contact_status.hasActiveContacts()) {
    robot.computeBaumgarteResidual(contact_status, dtau, contact_status.contactPoints(), data.C());
  }
  double l1norm_ref = dtau * data.IDC().lpNorm<1>();
  if (robot.hasFloatingBase()) {
     l1norm_ref += dtau * s.u_passive.lpNorm<1>();
  }
  double squarednorm_ref = dtau * dtau * data.IDC().squaredNorm();
  if (robot.hasFloatingBase()) {
     squarednorm_ref += dtau * dtau * s.u_passive.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


TEST_F(ContactDynamicsTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(fixed_base_urdf, contact_frames);
  contact_status.setContactStatus({false});
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondensing(robot, contact_status);
  testExpansionPrimal(robot, contact_status);
  testExpansionDual(robot, contact_status);
  testIntegration(robot, contact_status);
  testComputeResidual(robot, contact_status);
  contact_status.setContactStatus({true});
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondensing(robot, contact_status);
  testExpansionPrimal(robot, contact_status);
  testExpansionDual(robot, contact_status);
  testIntegration(robot, contact_status);
  testComputeResidual(robot, contact_status);
}


TEST_F(ContactDynamicsTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(floating_base_urdf, contact_frames);
  contact_status.setContactStatus({false, false, false, false});
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondensing(robot, contact_status);
  testExpansionPrimal(robot, contact_status);
  testExpansionDual(robot, contact_status);
  testIntegration(robot, contact_status);
  testComputeResidual(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  testLinearizeInverseDynamics(robot, contact_status);
  testLinearizeContactConstraints(robot, contact_status);
  testLinearizeContactDynamics(robot, contact_status);
  testCondensing(robot, contact_status);
  testExpansionPrimal(robot, contact_status);
  testExpansionDual(robot, contact_status);
  testIntegration(robot, contact_status);
  testComputeResidual(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}