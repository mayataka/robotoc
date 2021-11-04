#include <memory>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/impulse/impulse_split_ocp.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_dynamics.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class ImpulseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test_computeKKTResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status);
  static void test_computeKKTSystem(Robot& robot, 
                                    const ImpulseStatus& impulse_status);
  static void test_evalOCP(Robot& robot, const ImpulseStatus& impulse_status);
};


void ImpulseSplitOCPTest::test_computeKKTResidual(Robot& robot, 
                                                  const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  const double kkt_error = ocp.KKTError(kkt_residual);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost = cost->linearizeImpulseCost(robot, impulse_status, cost_data, t, s, kkt_residual_ref);
  constraints->linearizeConstraints(robot, constraints_data, s, kkt_residual_ref);
  impulse_cost += constraints_data.logBarrier();
  ImpulseStateEquation state_equation(robot);
  ImpulseStateEquation::linearizeStateEquation(robot, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamics id(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_residual_ref);
  kkt_residual_ref.kkt_error = ocp.KKTError(kkt_residual);
  const double kkt_error_ref = kkt_residual_ref.KKTError()
                                + id.KKTError() + constraints_data.KKTError();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_DOUBLE_EQ(impulse_cost, ocp.stageCost());
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void ImpulseSplitOCPTest::test_computeKKTSystem(Robot& robot, 
                                                const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ocp.computeKKTSystem(robot, impulse_status, t, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost = cost->quadratizeImpulseCost(robot, impulse_status, cost_data, t, s, kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, constraints_data, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(constraints_data, s, kkt_matrix_ref, kkt_residual_ref);
  impulse_cost += constraints_data.logBarrier();
  ImpulseStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  state_equation.correctLinearizedStateEquation(robot, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamics id(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_residual_ref );
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix_ref, kkt_residual_ref);

  kkt_residual.kkt_error = 0;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  ocp.expandPrimal(s, d);
  id.expandPrimal(d_ref);
  constraints->expandSlackAndDual(constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.expandDual(d_next, d);
  id.expandDual(d_next, d_ref);
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


void ImpulseSplitOCPTest::test_evalOCP(Robot& robot, const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  const auto d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ocp.evalOCP(robot, impulse_status, t, s, s_prev.q, s_prev.v, kkt_residual);
  const double impulse_cost = ocp.stageCost();
  const double constraint_violation = ocp.constraintViolation(kkt_residual);

  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost_ref = cost->evalImpulseCost(robot, impulse_status, cost_data, t, s);
  constraints->evalConstraint(robot, constraints_data, s);
  impulse_cost_ref +=  constraints_data.logBarrier();
  EXPECT_DOUBLE_EQ(impulse_cost, impulse_cost_ref);
  ImpulseStateEquation::evalStateEquation(robot, s, s_prev.q, s_prev.v, kkt_residual_ref);
  ImpulseDynamics id(robot);
  id.evalImpulseDynamics(robot, impulse_status, s);
  double constraint_violation_ref = 0;
  constraint_violation_ref += kkt_residual_ref.constraintViolation();
  constraint_violation_ref += constraints_data.constraintViolation();
  constraint_violation_ref += id.constraintViolation();
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}


TEST_F(ImpulseSplitOCPTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
}


TEST_F(ImpulseSplitOCPTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}