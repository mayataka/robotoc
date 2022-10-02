#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"
#include "robotoc/ocp/terminal_stage.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"
#include "contact_sequence_factory.hpp"


namespace robotoc {

class TerminalStageTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    grid_info.type = GridType::Terminal;
    grid_info.switching_constraint = false;
    grid_info.dt = 0.;
    N = 10;
    max_num_impulse = 10;
    t0 = 0.0;
    event_period = 0.1;
    t = grid_info.t;
  }

  virtual void TearDown() {
  }

  GridInfo grid_info;
  int N, max_num_impulse;
  double t, t0, event_period;
};


TEST_P(TerminalStageTest, evalOCP) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const SplitSolution s = SplitSolution::Random(robot);

  TerminalStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalOCP(robot, grid_info, s, data, kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);  
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  robot.updateKinematics(s.q, s.v);
  data_ref.performance_index.cost = cost->evalTerminalCost(robot, data_ref.cost_data, grid_info, s);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));
}


TEST_P(TerminalStageTest, evalKKT) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_prev = SplitSolution::Random(robot);

  TerminalStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTResidual kkt_residual(robot);  
  SplitKKTMatrix kkt_matrix(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalKKT(robot, grid_info, s_prev.q, s, data, kkt_matrix, kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);  
  SplitKKTMatrix kkt_matrix_ref(robot);  
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  robot.updateKinematics(s.q, s.v);
  data_ref.performance_index.cost = cost->quadratizeTerminalCost(robot, data_ref.cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  linearizeTerminalStateEquation(robot, s_prev.q, s, data_ref.state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  data_ref.performance_index.dual_feasibility = kkt_residual_ref.dualFeasibility();
  data_ref.performance_index.kkt_error = kkt_residual_ref.KKTError();
  correctLinearizeTerminalStateEquation(data_ref.state_equation_data, kkt_matrix_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));

  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  auto d_ocp = d;
  stage.expandPrimal(grid_info, data, d);
  stage.expandDual(grid_info, data, d);
  correctCostateDirection(data_ref.state_equation_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, TerminalStageTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(),
                    testhelper::CreateRobotManipulator(0.01),
                    testhelper::CreateQuadrupedalRobot(),
                    testhelper::CreateQuadrupedalRobot(0.01))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}