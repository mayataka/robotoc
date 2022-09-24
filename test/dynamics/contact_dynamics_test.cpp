#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"
#include "robotoc/dynamics/contact_dynamics.hpp"

#include "robot_factory.hpp"
#include "contact_status_factory.hpp"


namespace robotoc {

class ContactDynamicsTest : public ::testing::TestWithParam<std::pair<Robot, ContactStatus>> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dt;
};


TEST_P(ContactDynamicsTest, residual) {
  auto robot = GetParam().first;
  const auto contact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamicsData data(robot);
  evalContactDynamics(robot, contact_status, s, data);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data_ref.C());
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
}


TEST_P(ContactDynamicsTest, linearize) {
  auto robot = GetParam().first;
  const auto contact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  auto kkt_residual = SplitKKTResidual::Random(robot, contact_status);
  auto kkt_residual_ref = kkt_residual;
  ContactDynamicsData data(robot);
  linearizeContactDynamics(robot, contact_status, s, data, kkt_residual);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data_ref.C());
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data_ref.dCdq(), data_ref.dCdv(), data_ref.dCda());
  kkt_residual_ref.lq() += data_ref.dIDdq().transpose() * s.beta;
  kkt_residual_ref.lv() += data_ref.dIDdv().transpose() * s.beta;
  kkt_residual_ref.la   += data_ref.dIDda.transpose() * s.beta;
  kkt_residual_ref.lf() -= data_ref.dCda() * s.beta;
  Eigen::VectorXd lu_full = Eigen::VectorXd::Zero(robot.dimv());
  lu_full.head(robot.dim_passive()).setZero();
  lu_full.tail(robot.dimu()) = kkt_residual_ref.lu;
  lu_full -= s.beta;
  kkt_residual_ref.lq() += data_ref.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv() += data_ref.dCdv().transpose() * s.mu_stack(); 
  kkt_residual_ref.la   += data_ref.dCda().transpose() * s.mu_stack();
  lu_full.head(robot.dim_passive()) += s.nu_passive;
  data_ref.lu_passive     = lu_full.head(robot.dim_passive());
  kkt_residual_ref.lu = lu_full.tail(robot.dimu());
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
}


TEST_P(ContactDynamicsTest, condense) {
  auto robot = GetParam().first;
  const auto contact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  ContactDynamicsData data(robot);
  auto kkt_residual = SplitKKTResidual::Random(robot, contact_status);
  linearizeContactDynamics(robot, contact_status, s, data, kkt_residual);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, data_ref.ID_full());
  data_ref.ID().noalias() -= s.u;
  robot.computeBaumgarteResidual(contact_status, data_ref.C());
  robot.RNEADerivatives(s.q, s.v, s.a, data_ref.dIDdq(), data_ref.dIDdv(), data_ref.dIDda);
  robot.computeBaumgarteDerivatives(contact_status, data_ref.dCdq(), data_ref.dCdv(), data_ref.dCda());
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
  condenseContactDynamics(robot, contact_status, dt, data, kkt_matrix, kkt_residual);
  robot.computeMJtJinv(data_ref.dIDda, data_ref.dCda(), data_ref.MJtJinv());
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC()     = data_ref.MJtJinv() * data_ref.IDC();
  Eigen::MatrixXd Qaaff = Eigen::MatrixXd::Zero(dimv+dimf, dimv+dimf);
  Qaaff.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Qaa;
  Qaaff.bottomRightCorner(dimf, dimf) = kkt_matrix_ref.Qff();
  data_ref.Qafqv() = - Qaaff * data_ref.MJtJinv_dIDCdqv();
  data_ref.Qafqv().bottomLeftCorner(dimf, dimv) -= kkt_matrix_ref.Qqf().transpose();
  Eigen::MatrixXd IO_mat = Eigen::MatrixXd::Zero(dimv+dimf, dimv);
  IO_mat.topRows(dimv).setIdentity();
  data_ref.Qafu_full() = Qaaff * data_ref.MJtJinv() * IO_mat;
  data_ref.la() = kkt_residual_ref.la;
  data_ref.lf() = - kkt_residual_ref.lf();
  data_ref.laf() -= Qaaff * data_ref.MJtJinv() * data_ref.IDC();
  kkt_matrix_ref.Qxx -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafqv();
  kkt_matrix_ref.Qxx.topRows(dimv) += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_dIDCdqv().bottomRows(dimf);
  Eigen::MatrixXd Qxu_full = Eigen::MatrixXd::Zero(2*dimv, dimv);
  Qxu_full.rightCols(dimu) = kkt_matrix_ref.Qxu;
  Qxu_full -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qafu_full();
  Qxu_full.topRows(dimv) -= kkt_matrix_ref.Qqf() * data_ref.MJtJinv().bottomLeftCorner(dimf, dimv);
  data_ref.Qxu_passive   = Qxu_full.leftCols(dim_passive);
  kkt_matrix_ref.Qxu = Qxu_full.rightCols(dimu);
  const Eigen::MatrixXd Quu_full = IO_mat.transpose() * data_ref.MJtJinv() * data_ref.Qafu_full();
  data_ref.Quu_passive_topRight = Quu_full.topRightCorner(dim_passive, dimu);
  kkt_matrix_ref.Quu       += Quu_full.bottomRightCorner(dimu, dimu);
  kkt_residual_ref.lx -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.laf();
  kkt_residual_ref.lq() += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_IDC().tail(dimf);
  Eigen::VectorXd lu_full = Eigen::VectorXd::Zero(dimv);
  lu_full.head(dim_passive) = s.nu_passive - s.beta.head(dim_passive);
  lu_full.tail(dimu)        = kkt_residual_ref.lu;
  lu_full += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.laf();
  data_ref.lu_passive     = lu_full.head(dim_passive);
  kkt_residual_ref.lu = lu_full.tail(dimu);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv) = dt * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix_ref.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix_ref.Fxx  -= OOIO_mat * data_ref.MJtJinv_dIDCdqv();
  const Eigen::MatrixXd Fxu_full = OOIO_mat * data_ref.MJtJinv() * IO_mat;
  kkt_matrix_ref.Fvu = Fxu_full.bottomRows(dimv).rightCols(dimu);
  kkt_residual_ref.Fx -= (OOIO_mat * data_ref.MJtJinv() * data_ref.IDC());

  data_ref.ha() = kkt_matrix_ref.ha;
  data_ref.hf() = - kkt_matrix_ref.hf();
  kkt_residual_ref.h -= data_ref.MJtJinv_IDC().dot(data_ref.haf());
  kkt_matrix_ref.hx -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.haf();
  kkt_matrix_ref.hq() += (1.0/dt) * kkt_matrix_ref.Qqf() * data_ref.MJtJinv_IDC().tail(dimf);
  Eigen::VectorXd hu_full = Eigen::VectorXd::Zero(dimv);
  hu_full.tail(dimu)        = kkt_matrix_ref.hu;
  hu_full += IO_mat.transpose() * data_ref.MJtJinv() * data_ref.haf();
  kkt_matrix_ref.hu = hu_full.tail(dimu);

  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  EXPECT_TRUE(kkt_matrix.Quu.isApprox(kkt_matrix.Quu.transpose()));

  auto d = SplitDirection::Random(robot, contact_status);
  auto d_ref = d;
  expandContactDynamicsPrimal(data, d);
  IO_mat.topRows(dimv).setIdentity();
  Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
  du_full.tail(robot.dimu()) = d_ref.du;
  d_ref.daf() = - data_ref.MJtJinv() * (data_ref.dIDCdqv() * d.dx - IO_mat * du_full + data_ref.IDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));

  const auto d_next = SplitDirection::Random(robot);
  const double dts = Eigen::VectorXd::Random(1)[0];
  expandContactDynamicsDual(dt, dts, data, d_next, d);
  if (robot.hasFloatingBase()) {
    Eigen::VectorXd du_full = Eigen::VectorXd::Zero(dimv);
    du_full.tail(dimu) = d_ref.du; 
    d_ref.dnu_passive = - (data_ref.lu_passive + data_ref.Qxu_passive.transpose() * d_ref.dx 
                            + data_ref.Quu_passive_topRight * d_ref.du
                            + (IO_mat.transpose() * data_ref.MJtJinv() * OOIO_mat.transpose() * d_next.dlmdgmm).head(dim_passive));
  }
  d_ref.dbetamu() = - data_ref.MJtJinv() * (data_ref.Qafqv() * d_ref.dx 
                                        + data_ref.Qafu_full() * du_full 
                                        + OOIO_mat.transpose() * d_next.dlmdgmm
                                        + data_ref.laf() + dts * data_ref.haf());
  EXPECT_TRUE(d.isApprox(d_ref));
}


constexpr double dt = 0.01;

INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, ContactDynamicsTest, 
  ::testing::Values(std::make_pair(testhelper::CreateRobotManipulator(dt), 
                                   testhelper::CreateRobotManipulator(dt).createContactStatus()), 
                    std::make_pair(testhelper::CreateRobotManipulator(dt), 
                                   testhelper::CreateActiveContactStatus(testhelper::CreateRobotManipulator(dt), dt)),
                    std::make_pair(testhelper::CreateQuadrupedalRobot(dt), 
                                   testhelper::CreateQuadrupedalRobot(dt).createContactStatus()), 
                    std::make_pair(testhelper::CreateQuadrupedalRobot(dt), 
                                   testhelper::CreateActiveContactStatus(testhelper::CreateQuadrupedalRobot(dt), dt)))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}