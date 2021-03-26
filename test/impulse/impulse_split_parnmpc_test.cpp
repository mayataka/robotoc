#include <memory>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_parnmpc.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class ImpulseSplitParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void testLinearizeOCP(Robot& robot, 
                               const ImpulseStatus& impulse_status);
  static void testComputeKKTResidual(Robot& robot, 
                                     const ImpulseStatus& impulse_status);
  static void testCostAndConstraintViolation(Robot& robot, 
                                             const ImpulseStatus& impulse_status);
};


void ImpulseSplitParNMPCTest::testLinearizeOCP(Robot& robot,
                                               const ImpulseStatus& impulse_status) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const SplitSolution s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  parnmpc.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  parnmpc.linearizeOCP(robot, impulse_status, t, s_prev.q, s_prev.v, s, s_next, kkt_matrix, kkt_residual);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  cost->computeImpulseCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  cost->computeImpulseCostHessian(robot, cost_data, t, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeImpulseBackwardEuler(robot, s_prev.q, s_prev.v, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseImpulseBackwardEuler(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsBackwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref);
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  parnmpc.computeCondensedPrimalDirection(robot, kkt_matrix, s, d);
  id.computeCondensedPrimalDirection(robot,  kkt_matrix_ref, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(parnmpc.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(parnmpc.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  parnmpc.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, d);
  id.computeCondensedDualDirection(robot, d_ref);
  Eigen::VectorXd dlmd_ref = d_ref.dlmd();
  if (robot.hasFloatingBase()) {
    d_ref.dlmd().head(6) = kkt_matrix_ref.Fqq_inv.transpose() * dlmd_ref.head(6);
  }
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  parnmpc.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void ImpulseSplitParNMPCTest::testComputeKKTResidual(Robot& robot, 
                                                     const ImpulseStatus& impulse_status) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const SplitSolution s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  parnmpc.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  parnmpc.computeKKTResidual(robot, impulse_status, t, s_prev.q, s_prev.v, s, s_next, kkt_matrix, kkt_residual);
  const double kkt_error = parnmpc.squaredNormKKTResidual(kkt_residual);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  cost->computeImpulseCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual_ref);
  stateequation::linearizeImpulseBackwardEuler(robot, s_prev.q, s_prev.v, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsBackwardEuler id(robot);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref);
  const double kkt_error_ref = stateequation::squaredNormStateEuqationResidual(kkt_residual_ref)
                                + kkt_residual_ref.lx().squaredNorm()
                                + kkt_residual_ref.ldv.squaredNorm()
                                + kkt_residual_ref.lf().squaredNorm()
                                + id.squaredNormImpulseDynamicsResidual(kkt_residual_ref)
                                + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void ImpulseSplitParNMPCTest::testCostAndConstraintViolation(Robot& robot, 
                                                             const ImpulseStatus& impulse_status) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size = 0.3;
  parnmpc.initConstraints(robot, s);
  const double stage_cost = parnmpc.stageCost(robot, t, s, step_size);
  ImpulseSplitKKTResidual kkt_residual(robot);
  const double constraint_violation = parnmpc.constraintViolation(robot, impulse_status, t, s_prev.q, s_prev.v, s, kkt_residual);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeImpulseCost(robot, cost_data, t, s);
  stage_cost_ref += constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::computeImpulseBackwardEulerResidual(robot, s_prev.q, s_prev.v, 
                                                     s, kkt_residual_ref);
  ImpulseDynamicsBackwardEuler id(robot);
  id.computeImpulseDynamicsResidual(robot, impulse_status, s, kkt_residual_ref);
  double constraint_violation_ref = 0;
  constraint_violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += id.l1NormImpulseDynamicsResidual(kkt_residual_ref);
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}


TEST_F(ImpulseSplitParNMPCTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  testLinearizeOCP(robot, impulse_status);
  testComputeKKTResidual(robot, impulse_status);
  testCostAndConstraintViolation(robot, impulse_status);
  impulse_status.activateImpulse(0);
  testLinearizeOCP(robot, impulse_status);
  testComputeKKTResidual(robot, impulse_status);
  testCostAndConstraintViolation(robot, impulse_status);
}


TEST_F(ImpulseSplitParNMPCTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  testLinearizeOCP(robot, impulse_status);
  testComputeKKTResidual(robot, impulse_status);
  testCostAndConstraintViolation(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  testLinearizeOCP(robot, impulse_status);
  testComputeKKTResidual(robot, impulse_status);
  testCostAndConstraintViolation(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}