#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/impulse_state_equation.hpp"
#include "robotoc/dynamics/impulse_dynamics.hpp"
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
    grid_info.type = GridType::Impulse;
    grid_info.switching_constraint = false;
    grid_info.dt = 0.;
    grid_info.impulse_index = 0;

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


TEST_P(ImpactStageTest, evalOCP) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const auto impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index);
  const auto s = SplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);

  ImpactStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalOCP(robot, grid_info, s, s_next, data, kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_residual_ref.setContactDimension(impulse_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  data_ref.performance_index.cost = cost->evalImpulseCost(robot, impulse_status, data_ref.cost_data, grid_info, s);
  constraints->evalConstraint(robot, impulse_status, data_ref.constraints_data, s);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  evalImpulseStateEquation(robot, s, s_next, kkt_residual_ref);
  evalImpulseDynamics(robot, impulse_status, s, data_ref.contact_dynamics_data);
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));
}


TEST_P(ImpactStageTest, evalKKT) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const auto impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index);
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot, impulse_status);
  const auto s_next = SplitSolution::Random(robot);

  ImpactStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalKKT(robot, grid_info, s_prev.q, s, s_next, data, kkt_matrix, kkt_residual);

  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_matrix_ref.setContactDimension(impulse_status.dimf());
  kkt_residual_ref.setContactDimension(impulse_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  data_ref.performance_index.cost = cost->quadratizeImpulseCost(robot, impulse_status, data_ref.cost_data, grid_info, s, 
                                                                kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, impulse_status, data_ref.constraints_data, s, kkt_residual_ref);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  linearizeImpulseStateEquation(robot, s_prev.q, s, s_next, data_ref.state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  linearizeImpulseDynamics(robot, impulse_status, s, data_ref.contact_dynamics_data, kkt_residual_ref);
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  data_ref.performance_index.dual_feasibility 
      = data_ref.dualFeasibility<1>() + kkt_residual_ref.dualFeasibility<1>();
  data_ref.performance_index.kkt_error = data_ref.KKTError() + kkt_residual_ref.KKTError();
  constraints->condenseSlackAndDual(impulse_status, data_ref.constraints_data, 
                                    kkt_matrix_ref, kkt_residual_ref);
  condenseImpulseDynamics(robot, impulse_status, data_ref.contact_dynamics_data, 
                          kkt_matrix_ref, kkt_residual_ref);
  correctLinearizeImpulseStateEquation(robot, s, s_next, data_ref.state_equation_data, 
                                       kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));

  SplitDirection d = SplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  auto d_ocp = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  stage.expandPrimal(grid_info, data, d);
  stage.expandDual(grid_info, data, d_next, d);
  d_ref.setContactDimension(impulse_status.dimf());
  expandImpulseDynamicsPrimal(data_ref.contact_dynamics_data, d_ref);
  constraints->expandSlackAndDual(impulse_status, data_ref.constraints_data, d_ref);
  expandImpulseDynamicsDual(data_ref.contact_dynamics_data, d_next, d_ref);
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