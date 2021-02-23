#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/forward_switching_constraint.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {

class ForwardSwitchingConstraintTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    dt1 = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt2 = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testLinearizeSwitchingConstraint(Robot& robot) const;
  void testComputeSwitchingConstraintResidual(Robot& robot) const;

  double dt1, dt2;
  std::string fixed_base_urdf, floating_base_urdf;
};


void ForwardSwitchingConstraintTest::testLinearizeSwitchingConstraint(Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const SplitSolution s = SplitSolution::Random(robot, impulse_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom();
  kkt_residual.la.setRandom();
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  SplitStateConstraintJacobian jac(robot);
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  auto jac_ref = jac;
  robot.updateKinematics(s.q);
  ForwardSwitchingConstraint switching_constraint(robot);
  switching_constraint.linearizeSwitchingConstraint(robot, impulse_status, dt1, dt2,
                                                    s, kkt_matrix, kkt_residual, jac);
  jac_ref.setImpulseStatus(impulse_status);
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactResidual(impulse_status, impulse_status.contactPoints(), kkt_residual_ref.P());
  robot.computeContactDerivative(impulse_status, kkt_matrix_ref.Pq());
  if (robot.hasFloatingBase()) {
    robot.dIntegratedConfiguration(s.q, dq, jac_ref.dintegrate_dq);
    robot.dIntegratedVelocity(s.q, dq, jac_ref.dintegrate_dv);
    jac_ref.Phiq() = kkt_matrix_ref.Pq() * jac_ref.dintegrate_dq;
    jac_ref.Phiv() = (dt1+dt2) * kkt_matrix_ref.Pq() * jac_ref.dintegrate_dv;
    jac_ref.Phia() = (dt1*dt2) * kkt_matrix_ref.Pq() * jac_ref.dintegrate_dv;
  }
  else {
    jac_ref.Phiq() = kkt_matrix_ref.Pq();
    jac_ref.Phiv() = (dt1+dt2) * kkt_matrix_ref.Pq();
    jac_ref.Phia() = (dt1*dt2) * kkt_matrix_ref.Pq();
  }
  kkt_residual_ref.lq() += jac_ref.Phiq().transpose() * s.xi_stack();
  kkt_residual_ref.lv() += jac_ref.Phiv().transpose() * s.xi_stack();
  kkt_residual_ref.la   += jac_ref.Phia().transpose() * s.xi_stack();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(jac.isApprox(jac_ref));
  const double kkt = switching_constraint.squaredNormSwitchingConstraintResidual(kkt_residual);
  const double kkt_ref = kkt_residual_ref.P().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt, kkt_ref);
  const double l1 = switching_constraint.l1NormSwitchingConstraintResidual(kkt_residual);
  const double l1_ref = kkt_residual_ref.P().lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1, l1_ref);
}


void ForwardSwitchingConstraintTest::testComputeSwitchingConstraintResidual(Robot& robot) const {
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  const SplitSolution s = SplitSolution::Random(robot, impulse_status);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto kkt_residual_ref = kkt_residual;
  robot.updateKinematics(s.q);
  ForwardSwitchingConstraint switching_constraint(robot);
  switching_constraint.computeSwitchingConstraintResidual(robot, impulse_status, 
                                                          dt1, dt2, s, kkt_residual);
  const Eigen::VectorXd dq = (dt1+dt2) * s.v + (dt1*dt2) * s.a;
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.integrateConfiguration(s.q, dq, 1.0, q);
  robot.updateKinematics(q);
  robot.computeContactResidual(impulse_status, impulse_status.contactPoints(), kkt_residual_ref.P());
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  const double kkt = switching_constraint.squaredNormSwitchingConstraintResidual(kkt_residual);
  const double kkt_ref = kkt_residual_ref.P().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt, kkt_ref);
  const double l1 = switching_constraint.l1NormSwitchingConstraintResidual(kkt_residual);
  const double l1_ref = kkt_residual_ref.P().lpNorm<1>();
  EXPECT_DOUBLE_EQ(l1, l1_ref);
}


TEST_F(ForwardSwitchingConstraintTest, fixedbase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testLinearizeSwitchingConstraint(robot);
  testComputeSwitchingConstraintResidual(robot);
}


TEST_F(ForwardSwitchingConstraintTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testLinearizeSwitchingConstraint(robot);
  testComputeSwitchingConstraintResidual(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}