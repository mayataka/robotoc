#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/unconstr/parnmpc_terminal_stage.hpp"
#include "robotoc/dynamics/unconstr_dynamics.hpp"
#include "robotoc/dynamics/unconstr_state_equation.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class ParNMPCTerminalStageTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    robot = testhelper::CreateRobotManipulator();
    grid_info = GridInfo::Random();
    t = grid_info.t;
    dt = grid_info.dt;
    cost = testhelper::CreateCost(robot);
    constraints = testhelper::CreateConstraints(robot);
  }

  virtual void TearDown() {
  }

  Robot robot;
  GridInfo grid_info;
  double t, dt;
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
};



TEST_F(ParNMPCTerminalStageTest, evalOCP) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot);
  ParNMPCTerminalStage stage(robot, cost, constraints);
  UnconstrOCPData data = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  stage.evalOCP(robot, grid_info, s_prev.q, s_prev.v, s, data, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  PerformanceIndex performance_index_ref;
  UnconstrOCPData data_ref;
  data_ref.cost_data = cost->createCostFunctionData(robot);
  data_ref.constraints_data = constraints->createConstraintsData(robot, grid_info.stage);
  data_ref.unconstr_dynamics = UnconstrDynamics(robot);
  const auto contact_status = robot.createContactStatus();
  constraints->setSlackAndDual(robot, contact_status, grid_info, s, data_ref.constraints_data);
  robot.updateKinematics(s.q, s.v, s.a);
  performance_index_ref.cost = cost->evalStageCost(robot, contact_status, data_ref.cost_data, grid_info, s);
  performance_index_ref.cost += cost->evalTerminalCost(robot, data_ref.cost_data, grid_info, s);
  constraints->evalConstraint(robot, contact_status, grid_info, s, data_ref.constraints_data);
  performance_index_ref.cost_barrier = data_ref.constraints_data.logBarrier();
  evalUnconstrBackwardEuler(grid_info.dt, s_prev.q, s_prev.v, s, kkt_residual_ref);
  data_ref.unconstr_dynamics.evalUnconstrDynamics(robot, s);
  performance_index_ref.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  EXPECT_TRUE(data.performance_index.isApprox(performance_index_ref));
}


TEST_F(ParNMPCTerminalStageTest, evalKKT) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot);
  ParNMPCTerminalStage stage(robot, cost, constraints);
  UnconstrOCPData data = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  stage.evalKKT(robot, grid_info, s_prev.q, s_prev.v, s, data, kkt_matrix, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  PerformanceIndex performance_index_ref;
  UnconstrOCPData data_ref;
  data_ref.cost_data = cost->createCostFunctionData(robot);
  data_ref.constraints_data = constraints->createConstraintsData(robot, grid_info.stage);
  data_ref.unconstr_dynamics = UnconstrDynamics(robot);
  const auto contact_status = robot.createContactStatus();
  constraints->setSlackAndDual(robot, contact_status, grid_info, s, data_ref.constraints_data);
  robot.updateKinematics(s.q, s.v, s.a);
  performance_index_ref.cost = cost->quadratizeStageCost(robot, contact_status, data_ref.cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  performance_index_ref.cost += cost->quadratizeTerminalCost(robot, data_ref.cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, contact_status, grid_info, s, data_ref.constraints_data, kkt_residual_ref);
  performance_index_ref.cost_barrier = data_ref.constraints_data.logBarrier();
  linearizeUnconstrBackwardEulerTerminal(dt, s_prev.q, s_prev.v, s, kkt_matrix_ref, kkt_residual_ref);
  data_ref.unconstr_dynamics.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual_ref);
  performance_index_ref.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  performance_index_ref.dual_feasibility
      = data_ref.dualFeasibility<1>() + kkt_residual_ref.dualFeasibility<1>();
  performance_index_ref.kkt_error
      = data_ref.KKTError() + kkt_residual_ref.KKTError();
  constraints->condenseSlackAndDual(contact_status, grid_info, data_ref.constraints_data, kkt_matrix_ref, kkt_residual_ref);
  data_ref.unconstr_dynamics.condenseUnconstrDynamics(kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(performance_index_ref));

  SplitDirection d = SplitDirection::Random(robot);
  auto d_ref = d;
  stage.expandPrimalAndDual(grid_info, kkt_matrix, kkt_residual, data, d);
  constraints->expandSlackAndDual(contact_status, grid_info, d_ref, data_ref.constraints_data);
  data_ref.unconstr_dynamics.expandPrimal(d_ref);
  data_ref.unconstr_dynamics.expandDual(grid_info.dt, kkt_matrix_ref, kkt_residual_ref, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(stage.maxPrimalStepSize(data), constraints->maxSlackStepSize(data_ref.constraints_data));
  EXPECT_DOUBLE_EQ(stage.maxDualStepSize(data), constraints->maxDualStepSize(data_ref.constraints_data));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  stage.updatePrimal(robot, step_size, d, s_updated, data);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(data_ref.constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}