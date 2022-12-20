#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/impact_state_equation.hpp"
#include "robotoc/dynamics/impact_dynamics.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/ocp/impact_stage.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"
#include "contact_sequence_factory.hpp"



namespace robotoc {

class ImpactStageTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    grid_info.type = GridType::Impact;
    grid_info.switching_constraint = false;
    grid_info.dt = 0.;
    grid_info.impact_index = 0;

    N = 10;
    max_num_impact = 10;
    t0 = 0.0;
    event_period = 0.1;
    t = grid_info.t;
  }

  virtual void TearDown() {
  }

  GridInfo grid_info;
  int N, max_num_impact;
  double t, t0, event_period;
};


TEST_P(ImpactStageTest, evalOCP) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impact, t0, event_period);
  const auto impact_status = contact_sequence->impactStatus(grid_info.impact_index);
  const auto s = SplitSolution::Random(robot, impact_status);
  const auto s_next = SplitSolution::Random(robot);

  ImpactStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalOCP(robot, grid_info, s, s_next, data, kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_residual_ref.setContactDimension(impact_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  const Eigen::VectorXd v_after_impact = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impact);
  data_ref.performance_index.cost = cost->evalImpactCost(robot, impact_status, data_ref.cost_data, grid_info, s);
  constraints->evalConstraint(robot, impact_status, grid_info, s, data_ref.constraints_data);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  evalImpactStateEquation(robot, s, s_next, kkt_residual_ref);
  evalImpactDynamics(robot, impact_status, s, data_ref.contact_dynamics_data);
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));
}


TEST_P(ImpactStageTest, evalKKT) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impact, t0, event_period);
  const auto impact_status = contact_sequence->impactStatus(grid_info.impact_index);
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot, impact_status);
  const auto s_next = SplitSolution::Random(robot);

  ImpactStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalKKT(robot, grid_info, s_prev.q, s, s_next, data, kkt_matrix, kkt_residual);

  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_matrix_ref.setContactDimension(impact_status.dimf());
  kkt_residual_ref.setContactDimension(impact_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  const Eigen::VectorXd v_after_impact = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impact);
  data_ref.performance_index.cost = cost->quadratizeImpactCost(robot, impact_status, data_ref.cost_data, grid_info, s, 
                                                                kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, impact_status, grid_info, s, data_ref.constraints_data, kkt_residual_ref);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  linearizeImpactStateEquation(robot, s_prev.q, s, s_next, data_ref.state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  linearizeImpactDynamics(robot, impact_status, s, data_ref.contact_dynamics_data, kkt_residual_ref);
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  data_ref.performance_index.dual_feasibility 
      = data_ref.dualFeasibility<1>() + kkt_residual_ref.dualFeasibility<1>();
  data_ref.performance_index.kkt_error = data_ref.KKTError() + kkt_residual_ref.KKTError();
  constraints->condenseSlackAndDual(impact_status, grid_info, data_ref.constraints_data, 
                                    kkt_matrix_ref, kkt_residual_ref);
  condenseImpactDynamics(robot, impact_status, data_ref.contact_dynamics_data, 
                          kkt_matrix_ref, kkt_residual_ref);
  correctLinearizeImpactStateEquation(robot, s, s_next, data_ref.state_equation_data, 
                                       kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));

  SplitDirection d = SplitDirection::Random(robot, impact_status);
  auto d_ref = d;
  auto d_ocp = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  stage.expandPrimal(grid_info, data, d);
  stage.expandDual(grid_info, data, d_next, d);
  d_ref.setContactDimension(impact_status.dimf());
  expandImpactDynamicsPrimal(data_ref.contact_dynamics_data, d_ref);
  constraints->expandSlackAndDual(impact_status, grid_info, d_ref, data_ref.constraints_data);
  expandImpactDynamicsDual(data_ref.contact_dynamics_data, d_next, d_ref);
  correctCostateDirection(data_ref.state_equation_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, ImpactStageTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(0.01),
                    testhelper::CreateQuadrupedalRobot(0.01),
                    testhelper::CreateHumanoidRobot(0.01))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}