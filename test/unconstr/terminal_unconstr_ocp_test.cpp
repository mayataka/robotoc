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
#include "robotoc/unconstr/terminal_unconstr_ocp.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class TerminalUnconstrOCPTest : public ::testing::Test {
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


TEST_F(TerminalUnconstrOCPTest, evalOCP) {
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  TerminalUnconstrOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, grid_info, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.evalOCP(robot, grid_info, s, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  PerformanceIndex performance_index_ref;
  UnconstrOCPData data;
  data.cost_data = cost->createCostFunctionData(robot);
  data.constraints_data = constraints->createConstraintsData(robot, grid_info.time_stage);
  data.unconstr_dynamics = UnconstrDynamics(robot);
  const auto contact_status = robot.createContactStatus();
  // constraints->setSlackAndDual(robot, contact_status, data.constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  performance_index_ref.cost = cost->evalTerminalCost(robot, data.cost_data, grid_info, s);
  // constraints->evalConstraint(robot, contact_status, data.constraints_data, s);
  // performance_index_ref.cost_barrier = data.constraints_data.logBarrier();
  // performance_index_ref.primal_feasibility 
  //     = data.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  EXPECT_TRUE(ocp.getEval().isApprox(performance_index_ref));
}


TEST_F(TerminalUnconstrOCPTest, evalKKT) {
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  TerminalUnconstrOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, grid_info, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.evalKKT(robot, grid_info, s, kkt_matrix, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  PerformanceIndex performance_index_ref;
  UnconstrOCPData data;
  data.cost_data = cost->createCostFunctionData(robot);
  data.constraints_data = constraints->createConstraintsData(robot, grid_info.time_stage);
  data.unconstr_dynamics = UnconstrDynamics(robot);
  const auto contact_status = robot.createContactStatus();
  constraints->setSlackAndDual(robot, contact_status, data.constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  performance_index_ref.cost = cost->quadratizeTerminalCost(robot, data.cost_data, grid_info, s, kkt_residual_ref, kkt_matrix_ref);
  // constraints->linearizeConstraints(robot, contact_status, data.constraints_data, s, kkt_residual_ref);
  performance_index_ref.cost_barrier = data.constraints_data.logBarrier();
  unconstr::stateequation::linearizeForwardEulerTerminal(s, kkt_residual_ref);
  // data.unconstr_dynamics.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual_ref);
  // performance_index_ref.primal_feasibility 
  //     = data.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  performance_index_ref.dual_feasibility
      = data.dualFeasibility<1>() + kkt_residual_ref.dualFeasibility<1>();
  performance_index_ref.kkt_error
      = data.KKTError() + kkt_residual_ref.KKTError();
  // constraints->condenseSlackAndDual(contact_status, data.constraints_data, kkt_matrix_ref, kkt_residual_ref);
  // data.unconstr_dynamics.condenseUnconstrDynamics(kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(ocp.getEval().isApprox(performance_index_ref));

  SplitDirection d = SplitDirection::Random(robot);
  auto d_ref = d;
  ocp.expandPrimalAndDual(d);
  // constraints->expandSlackAndDual(contact_status, data.constraints_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(data.constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(data.constraints_data));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(data.constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}