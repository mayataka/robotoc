#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/contact_dynamics_data.hpp"
#include "robotoc/ocp/contact_dynamics.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ContactDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test_computeResidual(Robot& robot, const ContactStatus& contact_status) const;
  void test_linearize(Robot& robot, const ContactStatus& contact_status) const;
  void test_condense(Robot& robot, const ContactStatus& contact_status) const;

  double dt;
};


void ContactDynamicsTest::test_computeResidual(Robot& robot, const ContactStatus& contact_status) const {
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics cd(robot);
  cd.evalContactDynamics(robot, contact_status, s);
  const double l1norm = cd.constraintViolation();
  const double squarednorm = cd.KKTError();
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data.C());
  double l1norm_ref = data.IDC().lpNorm<1>();
  double squarednorm_ref = data.IDC().squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


void ContactDynamicsTest::test_linearize(Robot& robot, const ContactStatus& contact_status) const {
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics cd(robot);
  auto kkt_residual = SplitKKTResidual::Random(robot, contact_status);
  auto kkt_residual_ref = kkt_residual;
  cd.linearizeContactDynamics(robot, contact_status, s, kkt_residual);
  const double l1norm = cd.constraintViolation();
  const double squarednorm = cd.KKTError();
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data.C());
  robot.RNEADerivatives(s.q, s.v, s.a, data.dIDdq(), data.dIDdv(), data.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data.dCdq(), data.dCdv(), data.dCda());
  kkt_residual_ref.lq() += data.dIDdq().transpose() * s.beta;
  kkt_residual_ref.lv() += data.dIDdv().transpose() * s.beta;
  kkt_residual_ref.la   += data.dIDda.transpose() * s.beta;
  kkt_residual_ref.lf() -= data.dCda() * s.beta;
  Eigen::VectorXd lu_full = Eigen::VectorXd::Zero(robot.dimv());
  lu_full.head(robot.dim_passive()).setZero();
  lu_full.tail(robot.dimu()) = kkt_residual_ref.lu;
  lu_full -= s.beta;
  kkt_residual_ref.lq() += data.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += data.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.la   += data.dCda().transpose() * s.mu_stack();
  lu_full.head(robot.dim_passive()) += s.nu_passive;
  data.lu_passive     = lu_full.head(robot.dim_passive());
  kkt_residual_ref.lu = lu_full.tail(robot.dimu());
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  const double l1norm_ref = data.IDC().lpNorm<1>();
  const double squarednorm_ref = data.IDC().squaredNorm() + data.lu_passive.squaredNorm();
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  EXPECT_DOUBLE_EQ(squarednorm, squarednorm_ref);
}


void ContactDynamicsTest::test_condense(Robot& robot, const ContactStatus& contact_status) const {
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamics cd(robot);
  auto kkt_residual = SplitKKTResidual::Random(robot, contact_status);
  cd.linearizeContactDynamics(robot, contact_status, s, kkt_residual);
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data.ID_full());
  data.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data.C());
  robot.RNEADerivatives(s.q, s.v, s.a, data.dIDdq(), data.dIDdv(), data.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data.dCdq(), data.dCdv(), data.dCda());
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dim_passive = robot.dim_passive();
  const int dimf = contact_status.dimf();
  auto kkt_matrix = SplitKKTMatrix::Random(robot, contact_status);
  kkt_matrix.Qaa.setZero();
  kkt_matrix.Qaa.diagonal().setRandom();
  kkt_matrix.Fxx.setZero();
  kkt_matrix.Fvu.setZero();
  kkt_matrix.Fqq().setIdentity();
  if (robot.hasFloatingBase()) {
    kkt_matrix.Fqq().topLeftCorner(6, 6).setRandom();
  }
  kkt_matrix.Fqv() = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  cd.condenseContactDynamics(robot, contact_status, dt, kkt_matrix, kkt_residual);
  robot.computeMJtJinv(data.dIDda, data.dCda(), data.MJtJinv());
  data.MJtJinv_dIDCdqv() = data.MJtJinv() * data.dIDCdqv();
  data.MJtJinv_IDC()     = data.MJtJinv() * data.IDC();
  Eigen::MatrixXd Qaaff = Eigen::MatrixXd::Zero(dimv+dimf, dimv+dimf);
  Qaaff.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Qaa;
  Qaaff.bottomRightCorner(dimf, dimf) = kkt_matrix_ref.Qff();
  data.Qafqv() = - Qaaff * data.MJtJinv_dIDCdqv();
  data.Qafqv().bottomLeftCorner(dimf, dimv) -= kkt_matrix_ref.Qqf().transpose();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  data.Qafu_full() = Qaaff * data.MJtJinv() * IO_mat;
  data.la() = kkt_residual_ref.la;
  data.lf() = - kkt_residual_ref.lf();
  data.laf() -= Qaaff * data.MJtJinv() * data.IDC();
  kkt_matrix_ref.Qxx -= data.MJtJinv_dIDCdqv().transpose() * data.Qafqv();
  kkt_matrix_ref.Qxx.topRows(dimv) += kkt_matrix_ref.Qqf() * data.MJtJinv_dIDCdqv().bottomRows(dimf);
  Eigen::MatrixXd Qxu_full = Eigen::MatrixXd::Zero(2*dimv, dimv);
  Qxu_full.rightCols(dimu) = kkt_matrix_ref.Qxu;
  Qxu_full -= data.MJtJinv_dIDCdqv().transpose() * data.Qafu_full();
  Qxu_full.topRows(dimv) -= kkt_matrix_ref.Qqf() * data.MJtJinv().bottomLeftCorner(dimf, dimv);
  data.Qxu_passive   = Qxu_full.leftCols(dim_passive);
  kkt_matrix_ref.Qxu = Qxu_full.rightCols(dimu);
  const Eigen::MatrixXd Quu_full = IO_mat.transpose() * data.MJtJinv() * data.Qafu_full();
  data.Quu_passive_topRight = Quu_full.topRightCorner(dim_passive, dimu);
  kkt_matrix_ref.Quu       += Quu_full.bottomRightCorner(dimu, dimu);
  kkt_residual_ref.lx -= data.MJtJinv_dIDCdqv().transpose() * data.laf();
  kkt_residual_ref.lq() += kkt_matrix_ref.Qqf() * data.MJtJinv_IDC().tail(dimf);
  Eigen::VectorXd lu_full = Eigen::VectorXd::Zero(dimv);
  lu_full.head(dim_passive) = s.nu_passive - s.beta.head(dim_passive);
  lu_full.tail(dimu)        = kkt_residual_ref.lu;
  lu_full += IO_mat.transpose() * data.MJtJinv() * data.laf();
  data.lu_passive     = lu_full.head(dim_passive);
  kkt_residual_ref.lu = lu_full.tail(dimu);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix_ref.Fxx  -= OOIO_mat * data.MJtJinv_dIDCdqv();
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fvu = Fxu_full.bottomRows(dimv).rightCols(dimu);
  kkt_residual_ref.Fx -= (OOIO_mat * data.MJtJinv() * data.IDC());

  data.ha() = kkt_matrix_ref.ha;
  data.hf() = - kkt_matrix_ref.hf();
  kkt_residual_ref.h -= data.MJtJinv_IDC().dot(data.haf());
  kkt_matrix_ref.hx -= data.MJtJinv_dIDCdqv().transpose() * data.haf();
  kkt_matrix_ref.hq() += (1.0/dt) * kkt_matrix_ref.Qqf() * data.MJtJinv_IDC().tail(dimf);
  Eigen::VectorXd hu_full = Eigen::VectorXd::Zero(dimv);
  hu_full.tail(dimu)        = kkt_matrix_ref.hu;
  hu_full += IO_mat.transpose() * data.MJtJinv() * data.haf();
  kkt_matrix_ref.hu = hu_full.tail(dimu);

  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));

  auto d = SplitDirection::Random(robot, contact_status);
  auto d_ref = d;
  cd.expandPrimal(d);
  IO_mat.topRows(dimv).setIdentity();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  du_full.tail(robot.dimu()) = d_ref.du;
  d_ref.daf() = - data.MJtJinv() * (data.dIDCdqv() * d.dx - IO_mat * du_full + data.IDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));

  const auto d_next = SplitDirection::Random(robot);
  const double dts = Eigen::VectorXd::Random(1)[0];
  cd.expandDual(dt, dts, d_next, d);
  if (robot.hasFloatingBase()) {
    Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
    du_full.tail(dimu) = d_ref.du; 
    d_ref.dnu_passive = - (data.lu_passive + data.Qxu_passive.transpose() * d_ref.dx 
                            + data.Quu_passive_topRight * d_ref.du
                            + (IO_mat.transpose() * data.MJtJinv() * OOIO_mat.transpose() * d_next.dlmdgmm).head(dim_passive));
  }
  d_ref.dbetamu() = - data.MJtJinv() * (data.Qafqv() * d_ref.dx 
                                        + data.Qafu_full() * du_full 
                                        + OOIO_mat.transpose() * d_next.dlmdgmm
                                        + data.laf() + dts * data.haf());
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ContactDynamicsTest, robotManipulator) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random());
  }
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
  contact_status.activateContact(0);
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
}


TEST_F(ContactDynamicsTest, quadrupedalRobot) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random());
  }
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
}


TEST_F(ContactDynamicsTest, humanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.setContactPlacement(i, SE3::Random());
  }
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test_computeResidual(robot, contact_status);
  test_linearize(robot, contact_status);
  test_condense(robot, contact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}