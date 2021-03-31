#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class SplitParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testLinearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                        const bool switching_constraint=false) const;
  void testComputeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                              const bool switching_constraint=false) const;
  void testCostAndConstraintViolation(Robot& robot, const ContactStatus& contact_status, 
                                      const bool switching_constraint=false) const;

  double t, dt;
};


void SplitParNMPCTest::testLinearizeOCP(Robot& robot, 
                                        const ContactStatus& contact_status, 
                                        const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = SplitSolution::Random(robot, contact_status, impulse_status);
  }
  else {
    stmp = SplitSolution::Random(robot, contact_status);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  if (switching_constraint) {
    parnmpc.linearizeOCP(robot, contact_status, t, dt, s_prev.q, s_prev.v, s, 
                         s_next, kkt_matrix, kkt_residual, impulse_status);
  }
  else {
    parnmpc.linearizeOCP(robot, contact_status, t, dt, s_prev.q, s_prev.v, s, 
                         s_next, kkt_matrix, kkt_residual);
  }
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  cost->computeStageCostHessian(robot, cost_data, t, dt, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dt, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeBackwardEuler(robot, dt, s_prev.q, s_prev.v, s, 
                                        s_next, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseBackwardEuler(robot, dt, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, s, 
                                                      kkt_matrix_ref, kkt_residual_ref);
  }
  cd.condenseContactDynamics(robot, contact_status, dt, kkt_matrix_ref, kkt_residual_ref, false);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d;
  if (switching_constraint) {
    d = SplitDirection::Random(robot, contact_status, impulse_status);
  }
  else {
    d = SplitDirection::Random(robot, contact_status);
  }
  auto d_ref = d;
  parnmpc.computeCondensedPrimalDirection(robot, dt, s, d);
  cd.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(parnmpc.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(parnmpc.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  parnmpc.computeCondensedDualDirection(robot, dt, kkt_matrix, kkt_residual, d);
  cd.computeCondensedDualDirection(robot, dt, kkt_matrix, kkt_residual, d.dgmm(), d_ref);
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


void SplitParNMPCTest::testComputeKKTResidual(Robot& robot, 
                                              const ContactStatus& contact_status, 
                                              const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = SplitSolution::Random(robot, contact_status, impulse_status);
  }
  else {
    stmp = SplitSolution::Random(robot, contact_status);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  if (switching_constraint) {
    parnmpc.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s_prev.v, 
                               s, s_next, kkt_matrix, kkt_residual, impulse_status);
  }
  else {
    parnmpc.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s_prev.v, 
                               s, s_next, kkt_matrix, kkt_residual);
  }
  const double kkt_error = parnmpc.squaredNormKKTResidual(kkt_residual, dt);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  stateequation::linearizeBackwardEuler(robot, dt, s_prev.q, s_prev.v, s, s_next, 
                                        kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, s, 
                                                      kkt_matrix_ref, kkt_residual_ref);
  }
  double kkt_error_ref = kkt_residual_ref.Fx().squaredNorm()
                         + kkt_residual_ref.lx().squaredNorm()
                         + kkt_residual_ref.la.squaredNorm()
                         + kkt_residual_ref.lf().squaredNorm()
                         + kkt_residual_ref.lu().squaredNorm()
                         + cd.squaredNormContactDynamicsResidual(dt)
                         + dt * dt * constraints->squaredNormPrimalAndDualResidual(constraints_data);
  if (robot.hasFloatingBase()) {
    kkt_error_ref += kkt_residual_ref.lu_passive.squaredNorm();
  }
  if (switching_constraint) {
    kkt_error_ref += switchingconstraint::squaredNormSwitchingConstraintResidual(kkt_residual_ref);
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitParNMPCTest::testCostAndConstraintViolation(Robot& robot, 
                                                      const ContactStatus& contact_status, 
                                                      const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = SplitSolution::Random(robot, contact_status, impulse_status);
  }
  else {
    stmp = SplitSolution::Random(robot, contact_status);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitParNMPC parnmpc(robot, cost, constraints);
  const double step_size = 0.3;
  parnmpc.initConstraints(robot, 10, s);
  const double stage_cost = parnmpc.stageCost(robot, t, dt, s, step_size);
  SplitKKTResidual kkt_residual(robot);
  double constraint_violation;
  if (switching_constraint) {
    constraint_violation = parnmpc.constraintViolation(robot, contact_status, t, dt, 
                                                       s_prev.q, s_prev.v, s, kkt_residual, impulse_status);
  }
  else {
    constraint_violation = parnmpc.constraintViolation(robot, contact_status, t, dt, 
                                                       s_prev.q, s_prev.v, s, kkt_residual);
  }
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeStageCost(robot, cost_data, t, dt, s);
  stage_cost_ref += dt * constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::computeBackwardEulerResidual(robot, dt, s_prev.q, s_prev.v, s, kkt_residual_ref);
  ContactDynamics cd(robot);
  cd.computeContactDynamicsResidual(robot, contact_status, s);
  if (switching_constraint) {
    switchingconstraint::computeSwitchingConstraintResidual(robot, impulse_status, 
                                                            kkt_residual_ref);
  }
  double constraint_violation_ref = 0;
  constraint_violation_ref += dt * constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += cd.l1NormContactDynamicsResidual(dt);
  if (switching_constraint) {
    constraint_violation_ref += switchingconstraint::l1NormSwitchingConstraintResidual(kkt_residual_ref);
  }
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitParNMPCTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  testLinearizeOCP(robot, contact_status);
  testComputeKKTResidual(robot, contact_status);
  testCostAndConstraintViolation(robot, contact_status);
  contact_status.activateContact(0);
  testLinearizeOCP(robot, contact_status);
  testComputeKKTResidual(robot, contact_status);
  testCostAndConstraintViolation(robot, contact_status);
  testLinearizeOCP(robot, contact_status, true);
  testComputeKKTResidual(robot, contact_status, true);
  testCostAndConstraintViolation(robot, contact_status, true);
}


TEST_F(SplitParNMPCTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  testLinearizeOCP(robot, contact_status);
  testComputeKKTResidual(robot, contact_status);
  testCostAndConstraintViolation(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  testLinearizeOCP(robot, contact_status);
  testComputeKKTResidual(robot, contact_status);
  testCostAndConstraintViolation(robot, contact_status);
  testLinearizeOCP(robot, contact_status, true);
  testComputeKKTResidual(robot, contact_status, true);
  testCostAndConstraintViolation(robot, contact_status, true);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}