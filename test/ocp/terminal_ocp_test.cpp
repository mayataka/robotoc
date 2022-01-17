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

class TerminalOCPTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    t = grid_info.t;
  }

  virtual void TearDown() {
  }

  GridInfo grid_info;
  double t;
};


TEST_P(TerminalOCPTest, computeKKTResidual) {
  auto robot = GetParam();
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  SplitKKTMatrix kkt_matrix(robot);  
  ocp.computeKKTResidual(robot, grid_info, s_prev.q, s, kkt_matrix, kkt_residual);
  const double KKT = ocp.KKTError(kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTResidual kkt_residual_ref(robot);  
  SplitKKTMatrix kkt_matrix_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->linearizeTerminalCost(robot, cost_data, grid_info, s, kkt_residual_ref);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  kkt_residual_ref.kkt_error = ocp.KKTError(kkt_residual_ref);
  double KKT_ref = 0;
  KKT_ref += kkt_residual_ref.lx.squaredNorm();
  EXPECT_DOUBLE_EQ(KKT, KKT_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_P(TerminalOCPTest, computeKKTSystem) {
  auto robot = GetParam();
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const auto s = SplitSolution::Random(robot);
  const auto s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTMatrix kkt_matrix(robot);  
  SplitKKTResidual kkt_residual(robot);  
  ocp.computeKKTSystem(robot, grid_info, s_prev.q, s, kkt_matrix, kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTMatrix kkt_matrix_ref(robot);  
  SplitKKTResidual kkt_residual_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->quadratizeTerminalCost(robot, cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  // kkt_residual_ref.kkt_error = ocp.KKTError(kkt_residual_ref);
  kkt_residual.kkt_error = 0;
  state_equation.correctLinearizedStateEquation(kkt_matrix_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  ocp.expandDual(d);
  state_equation.correctCostateDirection(d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_P(TerminalOCPTest, evalOCP) {
  auto robot = GetParam();
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const auto s = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  ocp.evalOCP(robot, grid_info, s, kkt_residual);
  const double terminal_cost = ocp.terminalCost();
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost_ref = cost->evalTerminalCost(robot, cost_data, grid_info, s);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, TerminalOCPTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateQuadrupedalRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}