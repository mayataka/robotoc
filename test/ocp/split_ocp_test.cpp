#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/contact_dynamics.hpp"
#include "idocp/ocp/switching_constraint.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class SplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt_next = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testLinearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                        const bool switching_constraint=false) const;
  void testComputeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                              const bool switching_constraint=false) const;
  void testCostAndConstraintViolation(Robot& robot, const ContactStatus& contact_status, 
                                      const bool switching_constraint=false) const;

  double t, dt, dt_next;
};


void SplitOCPTest::testLinearizeOCP(Robot& robot, 
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
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitStateConstraintJacobian jac(robot);
  if (switching_constraint) {
    ocp.linearizeOCP(robot, contact_status, t, dt, s_prev.q, s, s_next, kkt_matrix, kkt_residual, 
                     impulse_status, dt_next, jac);
  }
  else {
    ocp.linearizeOCP(robot, contact_status, t, dt, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  }
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  SplitStateConstraintJacobian jac_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  cost->computeStageCostHessian(robot, cost_data, t, dt, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dt, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseForwardEuler(robot, dt, s, s_next.q, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    SwitchingConstraint sc(robot);
    sc.linearizeSwitchingConstraint(robot, impulse_status, dt, dt_next, s, 
                                    kkt_matrix_ref, kkt_residual_ref, jac_ref);
    cd.condenseContactDynamics(robot, contact_status, dt, kkt_matrix_ref, kkt_residual_ref);
    cd.condenseSwitchingConstraint(kkt_residual_ref, jac_ref);
    EXPECT_TRUE(jac.isApprox(jac_ref));
  }
  else {
    cd.condenseContactDynamics(robot, contact_status, dt, kkt_matrix_ref, kkt_residual_ref);
  }
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.P().isApprox(kkt_residual_ref.P()));
  SplitDirection d;
  if (switching_constraint) {
    d = SplitDirection::Random(robot, contact_status, impulse_status);
  }
  else {
    d = SplitDirection::Random(robot, contact_status);
  }
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  ocp.computeCondensedPrimalDirection(robot, dt, s, d);
  cd.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.computeCondensedDualDirection(robot, dt, kkt_matrix, kkt_residual, d_next, d);
  cd.computeCondensedDualDirection(robot, dt, kkt_matrix, kkt_residual, d_next.dgmm(), d_ref);
  Eigen::VectorXd dlmd_ref = d_ref.dlmd();
  if (robot.hasFloatingBase()) {
    d_ref.dlmd().head(6) = - kkt_matrix_ref.Fqq_prev_inv.transpose() * dlmd_ref.head(6);
  }
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void SplitOCPTest::testComputeKKTResidual(Robot& robot, 
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
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitStateConstraintJacobian jac(robot);
  if (switching_constraint) {
    ocp.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s, s_next, 
                           kkt_matrix, kkt_residual, 
                           impulse_status, dt_next, jac);
  }
  else {
    ocp.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s, s_next, 
                           kkt_matrix, kkt_residual);
  }
  const double kkt_error = ocp.squaredNormKKTResidual(kkt_residual, dt);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  SplitStateConstraintJacobian jac_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
    SwitchingConstraint sc(robot);
    sc.linearizeSwitchingConstraint(robot, impulse_status, dt, dt_next, s, 
                                    kkt_matrix_ref, kkt_residual_ref, jac_ref);
    EXPECT_TRUE(jac.isApprox(jac_ref));
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
    kkt_error_ref += kkt_residual_ref.P().squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitOCPTest::testCostAndConstraintViolation(Robot& robot, 
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
  SplitOCP ocp(robot, cost, constraints);
  const double step_size = 0.3;
  ocp.initConstraints(robot, 10, s);
  const double stage_cost = ocp.stageCost(robot, t, dt, s, step_size);
  SplitKKTResidual kkt_residual(robot);
  double constraint_violation;
  if (switching_constraint) {
    constraint_violation = ocp.constraintViolation(robot, contact_status, t, dt, 
                                                   s, s_next.q, s_next.v, 
                                                   kkt_residual, impulse_status, dt_next);
  }
  else {
    constraint_violation = ocp.constraintViolation(robot, contact_status, t, dt, 
                                                   s, s_next.q, s_next.v, kkt_residual);
  }
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeStageCost(robot, cost_data, t, dt, s);
  stage_cost_ref += dt * constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, s_next.q, 
                                             s_next.v, kkt_residual_ref);
  ContactDynamics cd(robot);
  cd.computeContactDynamicsResidual(robot, contact_status, s);
  if (switching_constraint) {
    kkt_residual_ref.setImpulseStatus(impulse_status);
    SwitchingConstraint sc(robot);
    sc.computeSwitchingConstraintResidual(robot, impulse_status, dt, dt_next, 
                                          s, kkt_residual_ref);
  }
  double constraint_violation_ref = 0;
  constraint_violation_ref += dt * constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += cd.l1NormContactDynamicsResidual(dt);
  if (switching_constraint) {
    constraint_violation_ref += kkt_residual_ref.P().lpNorm<1>();
  }
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitOCPTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
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


TEST_F(SplitOCPTest, floatingBase) {
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