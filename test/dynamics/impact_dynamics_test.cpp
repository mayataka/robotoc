#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"
#include "robotoc/dynamics/impact_dynamics.hpp"

#include "robot_factory.hpp"
#include "contact_status_factory.hpp"


namespace robotoc {

class ImpactDynamicsTest : public ::testing::TestWithParam<std::pair<Robot, ImpactStatus>> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }
};


TEST_P(ImpactDynamicsTest, residual) {
  auto robot = GetParam().first;
  const auto impact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ContactDynamicsData data(robot);
  evalImpactDynamics(robot, impact_status, s, data);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(impact_status.dimf());
  robot.setImpactForces(impact_status, s.f);
  robot.RNEAImpact(s.q, s.dv, data_ref.ID_full());
  robot.computeImpactVelocityResidual(impact_status, data_ref.C());
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
}


TEST_P(ImpactDynamicsTest, linearize) {
  auto robot = GetParam().first;
  const auto impact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ContactDynamicsData data(robot);
  auto kkt_residual = SplitKKTResidual::Random(robot, impact_status);
  auto kkt_residual_ref = kkt_residual;
  linearizeImpactDynamics(robot, impact_status, s, data, kkt_residual);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(impact_status.dimf());
  robot.setImpactForces(impact_status, s.f);
  robot.RNEAImpact(s.q, s.dv, data_ref.ID_full());
  robot.computeImpactVelocityResidual(impact_status, data_ref.C());
  robot.RNEAImpactDerivatives(s.q, s.dv, data_ref.dIDdq(), data_ref.dIDddv);
  robot.computeImpactVelocityDerivatives(impact_status, data_ref.dCdq(), data_ref.dCdv());
  kkt_residual_ref.lq().noalias() += data_ref.dIDdq().transpose() * s.beta;
  kkt_residual_ref.ldv.noalias()  += data_ref.dIDddv.transpose() * s.beta;
  kkt_residual_ref.lf().noalias() -= data_ref.dCdv() * s.beta;
  kkt_residual_ref.lq().noalias() += data_ref.dCdq().transpose() * s.mu_stack();
  kkt_residual_ref.lv().noalias() += data_ref.dCdv().transpose() * s.mu_stack();
  kkt_residual_ref.ldv.noalias()  += data_ref.dCdv().transpose() * s.mu_stack();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.IDC().isApprox(data_ref.IDC()));
}


TEST_P(ImpactDynamicsTest, condense) {
  auto robot = GetParam().first;
  const auto impact_status = GetParam().second;
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateKinematics(s.q, s.v+s.dv);
  ContactDynamicsData data(robot);
  auto kkt_residual = SplitKKTResidual::Random(robot, impact_status);
  linearizeImpactDynamics(robot, impact_status, s, data, kkt_residual);
  ContactDynamicsData data_ref(robot);
  data_ref.setContactDimension(impact_status.dimf());
  robot.setImpactForces(impact_status, s.f);
  robot.RNEAImpact(s.q, s.dv, data_ref.ID_full());
  robot.computeImpactVelocityResidual(impact_status, data_ref.C());
  robot.RNEAImpactDerivatives(s.q, s.dv, data_ref.dIDdq(), data_ref.dIDddv);
  robot.computeImpactVelocityDerivatives(impact_status, data_ref.dCdq(), data_ref.dCdv());
  const int dimv = robot.dimv();
  const int dimf = impact_status.dimf();
  // auto kkt_matrix = SplitKKTMatrix::Random(robot, impact_status);
  auto kkt_matrix = SplitKKTMatrix(robot);
  kkt_matrix.setContactDimension(impact_status.dimf());
  kkt_matrix.setRandom();
  kkt_matrix.Fxx.setZero();
  kkt_matrix.Qdvdv.setZero();
  kkt_matrix.Qdvdv.diagonal().setRandom();
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  condenseImpactDynamics(robot, impact_status, data, kkt_matrix, kkt_residual);
  robot.computeMJtJinv(data_ref.dIDddv, data_ref.dCdv(), data_ref.MJtJinv());
  data_ref.MJtJinv_dIDCdqv() = data_ref.MJtJinv() * data_ref.dIDCdqv();
  data_ref.MJtJinv_IDC()     = data_ref.MJtJinv() * data_ref.IDC();
  Eigen::MatrixXd Qdvdvff = Eigen::MatrixXd::Zero(dimv+dimf, dimv+dimf);
  Qdvdvff.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Qdvdv;
  Qdvdvff.bottomRightCorner(dimf, dimf) = kkt_matrix_ref.Qff();
  data_ref.Qdvfqv() = - Qdvdvff * data_ref.MJtJinv_dIDCdqv();
  data_ref.Qdvfqv().bottomLeftCorner(dimf, dimv) -= kkt_matrix_ref.Qqf().transpose();
  data_ref.ldvf().head(dimv) = kkt_residual_ref.ldv;
  data_ref.ldvf().tail(dimf) = - kkt_residual_ref.lf();
  data_ref.ldvf() -= Qdvdvff * data_ref.MJtJinv_IDC();
  kkt_matrix_ref.Qxx -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.Qdvfqv();
  kkt_matrix_ref.Qxx.topRows(dimv) += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_dIDCdqv().bottomRows(dimf);
  kkt_residual_ref.lx -= data_ref.MJtJinv_dIDCdqv().transpose() * data_ref.ldvf();
  kkt_residual_ref.lx.head(dimv) += kkt_matrix_ref.Qqf() * data_ref.MJtJinv_IDC().tail(dimf);
  Eigen::MatrixXd OOIO_mat = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  OOIO_mat.bottomLeftCorner(dimv, dimv).setIdentity();
  kkt_matrix_ref.Fvv().setIdentity();
  kkt_matrix_ref.Fxx -= OOIO_mat * data_ref.MJtJinv_dIDCdqv();
  kkt_residual_ref.Fx -= OOIO_mat * data_ref.MJtJinv_IDC();
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));
  EXPECT_TRUE(kkt_matrix.Qxx.isApprox(kkt_matrix.Qxx.transpose()));
  auto d = SplitDirection::Random(robot, impact_status);
  auto d_ref = d;
  expandImpactDynamicsPrimal(data, d);
  d_ref.ddvf() = - data_ref.MJtJinv() * (data_ref.dIDCdqv() * d_ref.dx + data_ref.IDC());
  d_ref.df().array() *= -1;
  EXPECT_TRUE(d.isApprox(d_ref));
  const auto d_next = SplitDirection::Random(robot);
  expandImpactDynamicsDual(data, d_next, d);
  d_ref.dbetamu() = - data_ref.MJtJinv() * (data_ref.Qdvfqv() * d.dx 
                                         + OOIO_mat.transpose() * d_next.dlmdgmm
                                         + data_ref.ldvf());
  EXPECT_TRUE(d.isApprox(d_ref));
}


constexpr double dt = 0.01;

INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, ImpactDynamicsTest, 
  ::testing::Values(std::make_pair(testhelper::CreateRobotManipulator(dt), 
                                   testhelper::CreateRobotManipulator(dt).createImpactStatus()), 
                    std::make_pair(testhelper::CreateRobotManipulator(dt), 
                                   testhelper::CreateActiveImpactStatus(testhelper::CreateRobotManipulator(dt), dt)),
                    std::make_pair(testhelper::CreateQuadrupedalRobot(dt), 
                                   testhelper::CreateQuadrupedalRobot(dt).createImpactStatus()), 
                    std::make_pair(testhelper::CreateQuadrupedalRobot(dt), 
                                   testhelper::CreateActiveImpactStatus(testhelper::CreateQuadrupedalRobot(dt), dt)))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}