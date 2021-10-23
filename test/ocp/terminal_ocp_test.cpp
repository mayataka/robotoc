#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class TerminalOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void test_computeKKTResidual(Robot& robot);
  static void test_computeKKTSystem(Robot& robot);
  static void test_evalOCP(Robot& robot);
};


void TerminalOCPTest::test_computeKKTResidual(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  SplitKKTMatrix kkt_matrix(robot);  
  ocp.computeKKTResidual(robot, t, s_prev.q, s, kkt_matrix, kkt_residual);
  const double KKT = ocp.KKTError(kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTResidual kkt_residual_ref(robot);  
  SplitKKTMatrix kkt_matrix_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->linearizeTerminalCost(robot, cost_data, t, s, kkt_residual_ref);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  double KKT_ref = 0;
  KKT_ref += kkt_residual_ref.lx.squaredNorm();
  EXPECT_DOUBLE_EQ(KKT, KKT_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void TerminalOCPTest::test_computeKKTSystem(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const auto s = SplitSolution::Random(robot);
  const auto s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTMatrix kkt_matrix(robot);  
  SplitKKTResidual kkt_residual(robot);  
  ocp.computeKKTSystem(robot, t, s_prev.q, s, kkt_matrix, kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTMatrix kkt_matrix_ref(robot);  
  SplitKKTResidual kkt_residual_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->quadratizeTerminalCost(robot, cost_data, t, s, kkt_residual_ref, kkt_matrix_ref);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquationAlongLieGroup(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  ocp.expandDual(d);
  state_equation.correctCostateDirection(d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void TerminalOCPTest::test_evalOCP(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const auto s = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  ocp.evalOCP(robot, t, s, kkt_residual);
  const double terminal_cost = ocp.terminalCost();
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost_ref = cost->computeTerminalCost(robot, cost_data, t, s);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
}


TEST_F(TerminalOCPTest, fixedBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_evalOCP(robot);
}


TEST_F(TerminalOCPTest, floatingBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_evalOCP(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}