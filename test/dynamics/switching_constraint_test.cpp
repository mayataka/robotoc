#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/dynamics/switching_constraint.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class SwitchingConstraintTest  : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    dt1 = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt2 = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = dt1;
  }

  virtual void TearDown() {
  }
  double dt1, dt2, dt;
};


TEST_P(SwitchingConstraintTest, eval) {
  auto robot = GetParam();
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const SplitSolution s = SplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q);
  SwitchingConstraintData data(robot);
  SplitKKTResidual kkt_residual(robot);
  auto kkt_residual_ref = kkt_residual;
  evalSwitchingConstraint(robot, impulse_status, data, dt1, dt2, s, kkt_residual);
  kkt_residual_ref.setSwitchingConstraintDimension(impulse_status.dimf());
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactPositionResidual(impulse_status, kkt_residual_ref.P());
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_P(SwitchingConstraintTest, linearize) {
  auto robot = GetParam();
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const SplitSolution s = SplitSolution::Random(robot, impulse_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx.setRandom();
  kkt_residual.la.setRandom();
  auto kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix(robot);
  auto kkt_matrix_ref = kkt_matrix;
  SwitchingConstraintData data(robot);
  auto data_ref = data;
  robot.updateKinematics(s.q);
  linearizeSwitchingConstraint(robot, impulse_status, data, dt1, dt2,
                               s, kkt_matrix, kkt_residual);
  data_ref.setDimension(impulse_status.dimf());
  kkt_matrix_ref.setSwitchingConstraintDimension(impulse_status.dimf());
  kkt_residual_ref.setSwitchingConstraintDimension(impulse_status.dimf());
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactPositionResidual(impulse_status, kkt_residual_ref.P());
  robot.computeContactPositionDerivative(impulse_status, data_ref.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, dq, data_ref.Pq(), kkt_matrix_ref.Phiq());
    robot.dIntegrateTransport_dv(s.q, dq, data_ref.Pq(), kkt_matrix_ref.Phiv());
    robot.dIntegrateTransport_dv(s.q, dq, data_ref.Pq(), kkt_matrix_ref.Phia());
    kkt_matrix_ref.Phiv().array() *= (dt1+dt2);
    kkt_matrix_ref.Phia().array() *= (dt1*dt2);
  }
  else {
    kkt_matrix_ref.Phiq() = data_ref.Pq();
    kkt_matrix_ref.Phiv() = (dt1+dt2) * data_ref.Pq();
    kkt_matrix_ref.Phia() = (dt1*dt2) * data_ref.Pq();
  }
  kkt_residual_ref.lx += kkt_matrix_ref.Phix().transpose() * s.xi_stack();
  kkt_residual_ref.la += kkt_matrix_ref.Phia().transpose() * s.xi_stack();
  const Eigen::VectorXd dqt = 2.0 * (s.v + dt1 * s.a);
  kkt_matrix_ref.Phit()  = data_ref.Pq() * dqt;
  kkt_residual_ref.h += s.xi_stack().dot(kkt_matrix_ref.Phit());
  const Eigen::VectorXd PqT_xi = data_ref.Pq().transpose() * s.xi_stack();
  kkt_matrix_ref.Qtt += 2.0 * PqT_xi.dot(s.a);
  kkt_matrix_ref.hv().noalias() += 2.0 * PqT_xi;
  kkt_matrix_ref.ha.noalias()  += (2.0*dt1) * PqT_xi;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, SwitchingConstraintTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(0.01),
                    testhelper::CreateQuadrupedalRobot(0.01))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}