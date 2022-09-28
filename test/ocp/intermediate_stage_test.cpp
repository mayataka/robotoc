#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/dynamics/state_equation.hpp"
#include "robotoc/dynamics/contact_dynamics.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"
#include "robotoc/dynamics/switching_constraint.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/intermediate_stage.hpp"
#include "robotoc/ocp/split_ocp.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"
#include "contact_sequence_factory.hpp"


namespace robotoc {

class IntermediateStageTest : public ::testing::TestWithParam<std::pair<Robot, bool>> {
protected:
  virtual void SetUp() {
    grid_info = GridInfo::Random();
    grid_info_next = GridInfo::Random();

    grid_info.type = GridType::Intermediate;
    grid_info.switching_constraint = false;
    grid_info.contact_phase = 0;
    grid_info.impulse_index = -1;
    grid_info.N_phase = 15;
    grid_info.dt_next = grid_info_next.dt;

    N = 10;
    max_num_impulse = 10;
    t0 = 0.0;
    event_period = 0.1;
    t = grid_info.t;
  }

  virtual void TearDown() {
  }

  GridInfo grid_info, grid_info_next;
  int N, max_num_impulse;
  double t, t0, event_period;
};


TEST_P(IntermediateStageTest, evalOCP) {
  auto robot = GetParam().first;
  const bool switching_constraint = GetParam().second;
  grid_info.switching_constraint = switching_constraint;
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const auto contact_status = contact_sequence->contactStatus(grid_info.contact_phase);
  auto s = SplitSolution::Random(robot, contact_status);
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    s.setSwitchingConstraintDimension(impulse_status.dimf());
    s.setRandom(robot);
  }
  const auto s_next = SplitSolution::Random(robot);

  IntermediateStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalOCP(robot, grid_info, s, s_next, data, kkt_residual);

  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_residual_ref.setContactDimension(contact_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  robot.updateKinematics(s.q, s.v, s.a);
  data_ref.performance_index.cost = cost->evalStageCost(robot, contact_status, data_ref.cost_data, grid_info, s);
  constraints->evalConstraint(robot, contact_status, data_ref.constraints_data, s);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  evalStateEquation(robot, grid_info.dt, s, s_next, kkt_residual_ref);
  evalContactDynamics(robot, contact_status, s, data_ref.contact_dynamics_data);
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    kkt_residual_ref.setSwitchingConstraintDimension(impulse_status.dimf());
    evalSwitchingConstraint(robot, impulse_status, data_ref.switching_constraint_data, grid_info.dt, grid_info.dt_next, s, kkt_residual_ref);
  }
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));

  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, contact_status, grid_info.time_stage, s);
  SplitKKTResidual kkt_residual_ocp(robot);
  SwitchingConstraintResidual switch_res(robot);
  double constraint_violation;
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    ocp.evalOCP(robot, contact_status, grid_info, s, s_next.q, s_next.v, 
                kkt_residual_ocp, impulse_status, grid_info_next, switch_res);
    constraint_violation = ocp.constraintViolation(kkt_residual_ocp, switch_res);
  }
  else {
    ocp.evalOCP(robot, contact_status, grid_info, s, s_next.q, s_next.v, kkt_residual_ocp);
    constraint_violation = ocp.constraintViolation(kkt_residual_ocp);
  }
  EXPECT_DOUBLE_EQ(ocp.stageCost(false), data_ref.performance_index.cost);
  EXPECT_DOUBLE_EQ(ocp.stageCost(true), data_ref.performance_index.cost+data_ref.performance_index.cost_barrier);
  EXPECT_DOUBLE_EQ(constraint_violation, data_ref.performance_index.primal_feasibility);
}


TEST_P(IntermediateStageTest, evalKKT) {
  auto robot = GetParam().first;
  const bool switching_constraint = GetParam().second;
  grid_info.switching_constraint = switching_constraint;
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  auto contact_sequence = testhelper::CreateContactSequence(robot, N, max_num_impulse, t0, event_period);
  const auto contact_status = contact_sequence->contactStatus(grid_info.contact_phase);
  auto s_prev = SplitSolution::Random(robot);
  auto s = SplitSolution::Random(robot, contact_status);
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    s.setSwitchingConstraintDimension(impulse_status.dimf());
    s.setRandom(robot);
  }
  const auto s_next = SplitSolution::Random(robot);

  IntermediateStage stage(cost, constraints, contact_sequence);
  auto data = stage.createData(robot);
  SplitKKTMatrix kkt_matrix(robot);  
  SplitKKTResidual kkt_residual(robot);  
  stage.initConstraints(robot, grid_info, s, data);
  stage.evalKKT(robot, grid_info, s_prev.q, s, s_next, data, kkt_matrix, kkt_residual);

  SplitKKTMatrix kkt_matrix_ref(robot);  
  SplitKKTResidual kkt_residual_ref(robot);  
  kkt_matrix_ref.setContactDimension(contact_status.dimf());
  kkt_residual_ref.setContactDimension(contact_status.dimf());
  auto data_ref = stage.createData(robot);
  stage.initConstraints(robot, grid_info, s, data_ref);

  robot.updateKinematics(s.q, s.v, s.a);
  data_ref.performance_index.cost = cost->quadratizeStageCost(robot, contact_status, data_ref.cost_data, grid_info, s, 
                                                              kkt_residual_ref, kkt_matrix_ref);
  kkt_residual_ref.h  = (1.0/grid_info.dt) * data_ref.performance_index.cost;
  kkt_matrix_ref.hx   = (1.0/grid_info.dt) * kkt_residual_ref.lx;
  kkt_matrix_ref.hu   = (1.0/grid_info.dt) * kkt_residual_ref.lu;
  kkt_matrix_ref.ha   = (1.0/grid_info.dt) * kkt_residual_ref.la;
  kkt_matrix_ref.hf() = (1.0/grid_info.dt) * kkt_residual_ref.lf();
  constraints->linearizeConstraints(robot, contact_status, data_ref.constraints_data, s, kkt_residual_ref);
  data_ref.performance_index.cost_barrier = data_ref.constraints_data.logBarrier();
  linearizeStateEquation(robot, grid_info.dt, s_prev.q, s, s_next, data_ref.state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  linearizeContactDynamics(robot, contact_status, s, data_ref.contact_dynamics_data, kkt_residual_ref);
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    linearizeSwitchingConstraint(robot, impulse_status, data_ref.switching_constraint_data, grid_info.dt, grid_info.dt_next, s, 
                                 kkt_matrix_ref, kkt_residual_ref);
  }
  data_ref.performance_index.primal_feasibility 
      = data_ref.primalFeasibility<1>() + kkt_residual_ref.primalFeasibility<1>();
  data_ref.performance_index.dual_feasibility 
      = data_ref.dualFeasibility<1>() + kkt_residual_ref.dualFeasibility<1>();
  data_ref.performance_index.kkt_error 
      = data_ref.KKTError() + kkt_residual_ref.KKTError();
  constraints->condenseSlackAndDual(contact_status, data_ref.constraints_data, 
                                    kkt_matrix_ref, kkt_residual_ref);
  condenseContactDynamics(robot, contact_status, grid_info.dt, 
                          data_ref.contact_dynamics_data, kkt_matrix_ref, kkt_residual_ref);
  correctLinearizeStateEquation(robot, grid_info.dt, s, s_next, 
                                data_ref.state_equation_data, kkt_matrix_ref, kkt_residual_ref);
  kkt_residual_ref.h        *= (1.0 / grid_info.N_phase);
  kkt_matrix_ref.hx.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix_ref.hu.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix_ref.fx.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix_ref.Qtt        *= 1.0 / (grid_info.N_phase * grid_info.N_phase);
  kkt_matrix_ref.Qtt_prev    = - kkt_matrix_ref.Qtt;
  if (switching_constraint) {
    kkt_matrix_ref.Phit().array() *= (1.0/grid_info.N_phase);
  }
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(data.performance_index.isApprox(data_ref.performance_index));

  SplitDirection d = SplitDirection::Random(robot, contact_status);
  d.setSwitchingConstraintDimension(kkt_residual.dims());
  d.dxi().setRandom();
  d.dts_next = Eigen::VectorXd::Random(1)[0];
  d.dts = Eigen::VectorXd::Random(1)[0];
  auto d_ref = d;
  auto d_ocp = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  const double dts = (d.dts_next - d.dts) / grid_info.N_phase;
  stage.expandPrimal(grid_info, data, d);
  stage.expandDual(grid_info, data, d_next, d);
  d_ref.setContactDimension(contact_status.dimf());
  expandContactDynamicsPrimal(data_ref.contact_dynamics_data, d_ref);
  constraints->expandSlackAndDual(contact_status, data_ref.constraints_data, d_ref);
  expandContactDynamicsDual(grid_info.dt, dts, data_ref.contact_dynamics_data, d_next, d_ref);
  correctCostateDirection(data_ref.state_equation_data, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));

  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, contact_status, grid_info.time_stage, s);
  SplitKKTMatrix kkt_matrix_ocp(robot);
  SplitKKTResidual kkt_residual_ocp(robot);
  SwitchingConstraintJacobian switch_jac_ocp(robot);
  SwitchingConstraintResidual switch_res_ocp(robot);
  double constraint_violation;
  if (switching_constraint) {
    const auto& impulse_status = contact_sequence->impulseStatus(grid_info.impulse_index+1);
    ocp.computeKKTSystem(robot, contact_status, grid_info, s_prev.q, s, s_next, kkt_matrix_ocp, kkt_residual_ocp, 
                         impulse_status, grid_info_next, switch_jac_ocp, switch_res_ocp);
    ocp.correctSTOSensitivities(kkt_matrix_ocp, kkt_residual_ocp, switch_jac_ocp, grid_info.N_phase);
    EXPECT_TRUE(kkt_matrix.Phiq().isApprox(switch_jac_ocp.Phiq()));
    EXPECT_TRUE(kkt_matrix.Phiv().isApprox(switch_jac_ocp.Phiv()));
    EXPECT_TRUE(kkt_matrix.Phiu().isApprox(switch_jac_ocp.Phiu()));
    EXPECT_TRUE(kkt_matrix.Phit().isApprox(switch_jac_ocp.Phit()));
    EXPECT_TRUE(kkt_residual.P().isApprox(switch_res_ocp.P()));
  }
  else {
    ocp.computeKKTSystem(robot, contact_status, grid_info, s_prev.q, s, s_next, kkt_matrix_ocp, kkt_residual_ocp);
    ocp.correctSTOSensitivities(kkt_matrix_ocp, kkt_residual_ocp, grid_info.N_phase);
  }
  kkt_matrix.setSwitchingConstraintDimension(0);
  kkt_residual.setSwitchingConstraintDimension(0);
  kkt_residual_ocp.kkt_error = 0.0;
  kkt_residual_ocp.cost = 0.0;
  kkt_residual_ocp.constraint_violation = 0.0;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ocp));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ocp));

  ocp.expandPrimal(contact_status, d_ocp);
  if (switching_constraint) {
    ocp.expandDual(grid_info, d_next, switch_jac_ocp, d_ocp, dts);
  }
  else {
    ocp.expandDual(grid_info, d_next, d_ocp, dts);
  }
  EXPECT_TRUE(d.isApprox(d_ocp));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, IntermediateStageTest, 
  ::testing::Values(std::make_pair(testhelper::CreateRobotManipulator(), false),
                    std::make_pair(testhelper::CreateRobotManipulator(0.01), false),
                    std::make_pair(testhelper::CreateRobotManipulator(0.01), true),
                    std::make_pair(testhelper::CreateQuadrupedalRobot(), false),
                    std::make_pair(testhelper::CreateQuadrupedalRobot(0.01), false),
                    std::make_pair(testhelper::CreateQuadrupedalRobot(0.01), true),
                    std::make_pair(testhelper::CreateHumanoidRobot(), false),
                    std::make_pair(testhelper::CreateHumanoidRobot(0.01), false),
                    std::make_pair(testhelper::CreateHumanoidRobot(0.01), true))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}