#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_switching_constraint_residual.hpp"
#include "idocp/ocp/split_switching_constraint_jacobian.hpp"
#include "idocp/ocp/switching_constraint.hpp"

#include "robot_factory.hpp"


namespace idocp {

class SwitchingConstraintTest : public ::testing::Test {
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

  void test_linearizeSwitchingConstraint(Robot& robot) const;
  void test_computeSwitchingConstraintResidual(Robot& robot) const;

  double dt1, dt2, dt;
};


void SwitchingConstraintTest::test_linearizeSwitchingConstraint(Robot& robot) const {
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
  SplitSwitchingConstraintJacobian jac(robot);
  SplitSwitchingConstraintResidual res(robot);
  auto jac_ref = jac;
  auto res_ref = res;
  robot.updateKinematics(s.q);
  switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, dt1, dt2,
                                                    s, kkt_matrix, kkt_residual, jac, res);
  jac_ref.setImpulseStatus(impulse_status);
  res_ref.setImpulseStatus(impulse_status);
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactPositionResidual(impulse_status, impulse_status.contactPoints(), res_ref.P());
  robot.computeContactPositionDerivative(impulse_status, jac_ref.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegrateTransport_dq(s.q, dq, jac_ref.Pq(), jac_ref.Phiq());
    robot.dIntegrateTransport_dv(s.q, dq, jac_ref.Pq(), jac_ref.Phiv());
    robot.dIntegrateTransport_dv(s.q, dq, jac_ref.Pq(), jac_ref.Phia());
    jac_ref.Phiv().array() *= (dt1+dt2);
    jac_ref.Phia().array() *= (dt1*dt2);
  }
  else {
    jac_ref.Phiq() = jac_ref.Pq();
    jac_ref.Phiv() = (dt1+dt2) * jac_ref.Pq();
    jac_ref.Phia() = (dt1*dt2) * jac_ref.Pq();
  }
  kkt_residual_ref.lx += jac_ref.Phix().transpose() * s.xi_stack();
  kkt_residual_ref.la += jac_ref.Phia().transpose() * s.xi_stack();
  res_ref.dq = s.v + dt1 * s.a;
  jac_ref.Phit().noalias() = jac_ref.Pq() * res_ref.dq;
  kkt_residual_ref.h = s.xi_stack().dot(jac_ref.Phit());
  kkt_matrix_ref.hq().setZero();
  kkt_matrix_ref.hv().noalias() = jac_ref.Pq().transpose() * s.xi_stack();
  kkt_matrix_ref.hu.setZero();
  res_ref.dq.noalias() = dt1 * jac_ref.Pq().transpose() * s.xi_stack();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(jac.isApprox(jac_ref));
  EXPECT_TRUE(res.isApprox(res_ref));
  EXPECT_TRUE(res.dq.isApprox(res_ref.dq));
  const double l2 = res.KKTError();
  const double l2_ref = res.P().squaredNorm();
  EXPECT_DOUBLE_EQ(l2, l2_ref);
  const double l1 = res.constraintViolation();
  const double l1_ref = res.P().lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1, l1_ref);
}


void SwitchingConstraintTest::test_computeSwitchingConstraintResidual(Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const SplitSolution s = SplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q);
  SplitSwitchingConstraintResidual res(robot);
  auto res_ref = res;
  switchingconstraint::computeSwitchingConstraintResidual(robot, impulse_status, 
                                                          dt1, dt2, s, res);

  res_ref.setImpulseStatus(impulse_status);
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactPositionResidual(impulse_status, impulse_status.contactPoints(), res_ref.P());
  EXPECT_TRUE(res.isApprox(res_ref));
  const double l2 = res.KKTError();
  const double l2_ref = res.P().squaredNorm();
  EXPECT_DOUBLE_EQ(l2, l2_ref);
  const double l1 = res.constraintViolation();
  const double l1_ref = res.P().lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1, l1_ref);
}


TEST_F(SwitchingConstraintTest, fixedbase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test_linearizeSwitchingConstraint(robot);
  test_computeSwitchingConstraintResidual(robot);
}


TEST_F(SwitchingConstraintTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test_linearizeSwitchingConstraint(robot);
  test_computeSwitchingConstraintResidual(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}