#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
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

class TerminalOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void testLinearizeOCP(Robot& robot);
  static void testTerminalCost(Robot& robot);
  static void testComputeKKTResidual(Robot& robot);
};


void TerminalOCPTest::testLinearizeOCP(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTMatrix kkt_matrix(robot);  
  SplitKKTResidual kkt_residual(robot);  
  ocp.linearizeOCP(robot, t, s_prev.q, s, kkt_matrix, kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTMatrix kkt_matrix_ref(robot);  
  SplitKKTResidual kkt_residual_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->quadratizeTerminalCost(robot, cost_data, t, s, kkt_residual_ref, kkt_matrix_ref);
  stateequation::linearizeForwardEulerTerminal(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseForwardEulerTerminal(robot, kkt_matrix_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  ocp.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, d);
  Eigen::VectorXd dlmd_ref = d_ref.dlmd();
  if (robot.hasFloatingBase()) {
    d_ref.dlmd().head(6) = - kkt_matrix_ref.Fqq_prev_inv.transpose() * dlmd_ref.head(6);
  }
  EXPECT_TRUE(d.isApprox(d_ref));
}


void TerminalOCPTest::testTerminalCost(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  const double terminal_cost = ocp.terminalCost(robot, t, s);
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost_ref = cost->computeTerminalCost(robot, cost_data, t, s);
  EXPECT_DOUBLE_EQ(terminal_cost, terminal_cost_ref);
}


void TerminalOCPTest::testComputeKKTResidual(Robot& robot) {
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_prev = SplitSolution::Random(robot);
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  TerminalOCP ocp(robot, cost, constraints);
  SplitKKTResidual kkt_residual(robot);  
  SplitKKTMatrix kkt_matrix(robot);  
  ocp.computeKKTResidual(robot, t, s_prev.q, s, kkt_matrix, kkt_residual);
  const double KKT = ocp.squaredNormKKTResidual(kkt_residual);
  robot.updateKinematics(s.q, s.v);
  SplitKKTResidual kkt_residual_ref(robot);  
  SplitKKTMatrix kkt_matrix_ref(robot);  
  robot.updateKinematics(s.q, s.v);
  auto cost_data = cost->createCostFunctionData(robot);
  const double terminal_cost = cost->linearizeTerminalCost(robot, cost_data, t, s, kkt_residual_ref);
  stateequation::linearizeForwardEulerTerminal(robot, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  double KKT_ref = 0;
  KKT_ref += kkt_residual_ref.lx.squaredNorm();
  EXPECT_DOUBLE_EQ(KKT, KKT_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(TerminalOCPTest, fixedBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testLinearizeOCP(robot);
  testTerminalCost(robot);
  testComputeKKTResidual(robot);
}


TEST_F(TerminalOCPTest, floatingBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testLinearizeOCP(robot);
  testTerminalCost(robot);
  testComputeKKTResidual(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}