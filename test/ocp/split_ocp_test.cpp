#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_ocp.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/state_equation.hpp"
#include "robotoc/ocp/contact_dynamics.hpp"
#include "robotoc/ocp/switching_constraint_residual.hpp"
#include "robotoc/ocp/switching_constraint_jacobian.hpp"
#include "robotoc/ocp/switching_constraint.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

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

  void test_computeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                               const bool switching_constraint=false) const;
  void test_computeKKTSystem(Robot& robot, const ContactStatus& contact_status, 
                             const bool switching_constraint=false) const;
  void test_evalOCP(Robot& robot, const ContactStatus& contact_status, 
                    const bool switching_constraint=false) const;

  double t, dt, dt_next;
};


void SplitOCPTest::test_computeKKTResidual(Robot& robot, 
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
  const auto s_prev = SplitSolution::Random(robot, contact_status);
  const auto s = stmp;
  const auto s_next = SplitSolution::Random(robot, contact_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SwitchingConstraintJacobian switch_jac(robot);
  SwitchingConstraintResidual switch_res(robot);
  if (switching_constraint) {
    ocp.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s, s_next, 
                           kkt_matrix, kkt_residual, 
                           impulse_status, dt_next, switch_jac, switch_res);
  }
  else {
    ocp.computeKKTResidual(robot, contact_status, t, dt, s_prev.q, s, s_next, 
                           kkt_matrix, kkt_residual);
  }
  const double kkt_error = ocp.KKTError(kkt_residual, dt);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost = cost->linearizeStageCost(robot, cost_data, t, dt, s, kkt_residual_ref);
  constraints->linearizeConstraints(robot, constraints_data, dt, s, kkt_residual_ref);
  stage_cost += dt * constraints_data.logBarrier();
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    SwitchingConstraintJacobian switch_jac_ref(robot);
    SwitchingConstraintResidual switch_res_ref(robot);
    switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, dt, dt_next, s, 
                                                      kkt_matrix_ref, kkt_residual_ref, 
                                                      switch_jac_ref, switch_res_ref);
    EXPECT_TRUE(switch_jac.isApprox(switch_jac_ref));
    EXPECT_TRUE(switch_res.isApprox(switch_res_ref));
  }
  const double kkt_error_ref = kkt_residual_ref.KKTError()
                                + (dt*dt) * cd.KKTError()
                                + (dt*dt) * constraints_data.KKTError();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_DOUBLE_EQ(stage_cost, ocp.stageCost());
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitOCPTest::test_computeKKTSystem(Robot& robot, 
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
  const auto s_prev = SplitSolution::Random(robot, contact_status);
  const auto s = stmp;
  const auto s_next = SplitSolution::Random(robot, contact_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SwitchingConstraintJacobian switch_jac(robot);
  SwitchingConstraintResidual switch_res(robot);
  if (switching_constraint) {
    ocp.computeKKTSystem(robot, contact_status, t, dt, s_prev.q, s, s_next, kkt_matrix, kkt_residual, 
                     impulse_status, dt_next, switch_jac, switch_res);
  }
  else {
    ocp.computeKKTSystem(robot, contact_status, t, dt, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  }
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  SwitchingConstraintJacobian switch_jac_ref(robot);
  SwitchingConstraintResidual switch_res_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost = cost->quadratizeStageCost(robot, cost_data, t, dt, s, kkt_residual_ref, kkt_matrix_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dt, s, kkt_matrix_ref, kkt_residual_ref);
  stage_cost += dt * constraints_data.logBarrier();
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquationAlongLieGroup(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dt, s, kkt_residual_ref);
  if (switching_constraint) {
    SwitchingConstraintJacobian switch_jac_ref(robot);
    SwitchingConstraintResidual switch_res_ref(robot);
    switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, dt, dt_next, s, 
                                                      kkt_matrix_ref, kkt_residual_ref, 
                                                      switch_jac_ref, switch_res_ref);
    cd.condenseContactDynamics(robot, contact_status, dt, s, kkt_matrix_ref, kkt_residual_ref);
    cd.condenseSwitchingConstraint(switch_jac_ref, switch_res_ref, kkt_matrix_ref);
    EXPECT_TRUE(switch_jac.isApprox(switch_jac_ref));
    EXPECT_TRUE(switch_res.isApprox(switch_res_ref));
  }
  else {
    cd.condenseContactDynamics(robot, contact_status, dt, s, kkt_matrix_ref, kkt_residual_ref);
  }
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
  const auto d_next = SplitDirection::Random(robot);
  const bool sto = false;
  ocp.expandPrimal(dt, s, d, sto);
  cd.expandPrimal(d_ref);
  constraints->expandSlackAndDual(constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.expandDual(dt, d_next, d, sto);
  cd.expandDual(dt, d_next, d_ref);
  state_equation.correctCostateDirection(d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void SplitOCPTest::test_evalOCP(Robot& robot, const ContactStatus& contact_status, 
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
  const auto s_prev = SplitSolution::Random(robot, contact_status);
  const auto s = stmp;
  const auto s_next = SplitSolution::Random(robot, contact_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  SplitOCP ocp(robot, cost, constraints);
  const double step_size = 0.3;
  ocp.initConstraints(robot, 10, s);
  SplitKKTResidual kkt_residual(robot);
  SwitchingConstraintResidual switch_res(robot);
  double constraint_violation;
  if (switching_constraint) {
    ocp.evalOCP(robot, contact_status, t, dt, s, s_next.q, s_next.v, 
                    kkt_residual, impulse_status, dt_next, switch_res);
    constraint_violation = ocp.constraintViolation(kkt_residual, dt, switch_res);
  }
  else {
    ocp.evalOCP(robot, contact_status, t, dt, s, s_next.q, s_next.v, 
                    kkt_residual);
    constraint_violation = ocp.constraintViolation(kkt_residual, dt);
  }
  const double stage_cost = ocp.stageCost();

  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = cost->evalStageCost(robot, cost_data, t, dt, s);
  constraints->evalConstraint(robot, constraints_data, s);
  stage_cost_ref += dt * constraints_data.logBarrier();
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  StateEquation state_equation(robot);
  state_equation.evalStateEquation(robot, dt, s, s_next.q, s_next.v, kkt_residual_ref);
  ContactDynamics cd(robot);
  cd.evalContactDynamics(robot, contact_status, s);
  double switch_violation_ref = 0;
  if (switching_constraint) {
    SwitchingConstraintResidual switch_res_ref(robot);
    switchingconstraint::evalSwitchingConstraint(robot, impulse_status, dt,  
                                                 dt_next, s, switch_res_ref);
    EXPECT_TRUE(switch_res.isApprox(switch_res_ref));
    switch_violation_ref = switch_res_ref.P().lpNorm<1>();
  }
  double constraint_violation_ref = 0;
  constraint_violation_ref += kkt_residual_ref.constraintViolation();
  constraint_violation_ref += dt * constraints_data.constraintViolation();
  constraint_violation_ref += dt * cd.constraintViolation();
  constraint_violation_ref += switch_violation_ref;
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitOCPTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  test_computeKKTResidual(robot, contact_status);
  test_computeKKTSystem(robot, contact_status);
  test_evalOCP(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test_computeKKTResidual(robot, contact_status);
  test_computeKKTSystem(robot, contact_status);
  test_evalOCP(robot, contact_status);
  test_computeKKTResidual(robot, contact_status, true);
  test_computeKKTSystem(robot, contact_status, true);
  test_evalOCP(robot, contact_status, true);
}


TEST_F(SplitOCPTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto contact_status = robot.createContactStatus();
  test_computeKKTResidual(robot, contact_status);
  test_computeKKTSystem(robot, contact_status);
  test_evalOCP(robot, contact_status);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  test_computeKKTResidual(robot, contact_status);
  test_computeKKTSystem(robot, contact_status);
  test_evalOCP(robot, contact_status);
  test_computeKKTResidual(robot, contact_status, true);
  test_computeKKTSystem(robot, contact_status, true);
  test_evalOCP(robot, contact_status, true);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}