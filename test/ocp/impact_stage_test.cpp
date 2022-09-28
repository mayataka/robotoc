#include <memory>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/ocp/impulse_split_ocp.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/impulse_state_equation.hpp"
#include "robotoc/dynamics/impulse_dynamics.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class ImpactStageTest : public ::testing::Test {
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


void ImpactStageTest::test_computeKKTResidual(Robot& robot, 
                                                  const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpactStage ocp(robot, cost, constraints);
  const auto grid_info = GridInfo::Random();
  const double t = grid_info.t;
  ocp.initConstraints(robot, impulse_status, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.computeKKTResidual(robot, impulse_status, grid_info, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  const double kkt_error = ocp.KKTError(kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactDimension(impulse_status.dimf());
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactDimension(impulse_status.dimf());
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, impulse_status, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost = cost->linearizeImpulseCost(robot, impulse_status, cost_data, grid_info, s, kkt_residual_ref);
  constraints->linearizeConstraints(robot, impulse_status, constraints_data, s, kkt_residual_ref);
  impulse_cost += constraints_data.logBarrier();
  StateEquationData state_equation_data(robot);
  linearizeImpulseStateEquation(robot, s_prev.q, s, s_next, state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamicsData id_data(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  linearizeImpulseDynamics(robot, impulse_status, s, id_data, kkt_residual_ref);
  kkt_residual_ref.kkt_error = ocp.KKTError(kkt_residual);
  const double kkt_error_ref = kkt_residual_ref.KKTError()
                                + id_data.KKTError() + constraints_data.KKTError();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_DOUBLE_EQ(impulse_cost, ocp.stageCost());
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void ImpactStageTest::test_computeKKTSystem(Robot& robot, 
                                                const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpactStage ocp(robot, cost, constraints);
  const auto grid_info = GridInfo::Random();
  const double t = grid_info.t;
  ocp.initConstraints(robot, impulse_status, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.computeKKTSystem(robot, impulse_status, grid_info, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactDimension(impulse_status.dimf());
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactDimension(impulse_status.dimf());
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, impulse_status, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost = cost->quadratizeImpulseCost(robot, impulse_status, cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, impulse_status, constraints_data, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(impulse_status, constraints_data, kkt_matrix_ref, kkt_residual_ref);
  impulse_cost += constraints_data.logBarrier();
  StateEquationData state_equation_data(robot);
  linearizeImpulseStateEquation(robot, s_prev.q, s, s_next, state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  correctLinearizeImpulseStateEquation(robot, s, s_next, state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamicsData id_data(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  linearizeImpulseDynamics(robot, impulse_status, s, id_data, kkt_residual_ref );
  condenseImpulseDynamics(robot, impulse_status, id_data, kkt_matrix_ref, kkt_residual_ref);

  kkt_residual.kkt_error = 0;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  ocp.expandPrimal(impulse_status, d);
  expandImpulseDynamicsPrimal(id_data, d_ref);
  constraints->expandSlackAndDual(impulse_status, constraints_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.expandDual(d_next, d);
  expandImpulseDynamicsDual(id_data, d_next, d_ref);
  correctCostateDirection(state_equation_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d, true);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void ImpactStageTest::test_evalOCP(Robot& robot, const ImpulseStatus& impulse_status) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot, impulse_status);
  const auto d = SplitDirection::Random(robot, impulse_status);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ImpactStage ocp(robot, cost, constraints);
  const auto grid_info = GridInfo::Random();
  const double t = grid_info.t;
  ocp.initConstraints(robot, impulse_status, s);
  SplitKKTResidual kkt_residual(robot);
  ocp.evalOCP(robot, impulse_status, grid_info, s, s_prev.q, s_prev.v, kkt_residual);
  const double impulse_cost = ocp.stageCost();
  const double constraint_violation = ocp.constraintViolation(kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactDimension(impulse_status.dimf());
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, -1);
  constraints->setSlackAndDual(robot, impulse_status, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double impulse_cost_ref = cost->evalImpulseCost(robot, impulse_status, cost_data, grid_info, s);
  constraints->evalConstraint(robot, impulse_status, constraints_data, s);
  impulse_cost_ref +=  constraints_data.logBarrier();
  EXPECT_DOUBLE_EQ(impulse_cost, impulse_cost_ref);
  StateEquationData state_equation_data(robot);
  evalImpulseStateEquation(robot, s, s_prev.q, s_prev.v, kkt_residual_ref);
  ContactDynamicsData id_data(robot);
  evalImpulseDynamics(robot, impulse_status, s, id_data);
  double constraint_violation_ref = 0;
  constraint_violation_ref += kkt_residual_ref.constraintViolation();
  constraint_violation_ref += constraints_data.constraintViolation();
  constraint_violation_ref += id_data.constraintViolation();
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}


TEST_F(ImpactStageTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test_computeKKTSystem(robot, impulse_status);
  test_computeKKTResidual(robot, impulse_status);
  test_evalOCP(robot, impulse_status);
}


TEST_F(ImpactStageTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
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


TEST_F(ImpactStageTest, humanoidRobot) {
  const double dt = 0.001;
  auto robot = testhelper::CreateHumanoidRobot(dt);
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