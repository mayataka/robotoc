#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/unconstr/split_unconstr_parnmpc.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class SplitUnconstrParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    cost = testhelper::CreateCost(robot);
    constraints = testhelper::CreateConstraints(robot);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
};


TEST_F(SplitUnconstrParNMPCTest, computeKKTSystem) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitUnconstrParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  parnmpc.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  parnmpc.computeKKTSystem(robot, t, dt, s_prev.q, s_prev.v, s, s_next, kkt_matrix, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost = cost->quadratizeStageCost(robot, cost_data, t, dt, s, kkt_residual_ref, kkt_matrix_ref);
  constraints->linearizeConstraints(robot, constraints_data, dt, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(constraints_data, dt, s, kkt_matrix_ref, kkt_residual_ref);
  stage_cost += dt * constraints_data.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(dt, s_prev.q, s_prev.v, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  UnconstrDynamics ud(robot);
  ud.linearizeUnconstrDynamics(robot, dt, s, kkt_residual_ref);
  ud.condenseUnconstrDynamics(kkt_matrix_ref, kkt_residual_ref);

  kkt_residual.kkt_error = 0;
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot);
  auto d_ref = d;
  parnmpc.expandPrimalAndDual(dt, s, kkt_matrix, kkt_residual, d);
  constraints->expandSlackAndDual(constraints_data, s, d_ref);
  ud.expandPrimal(d_ref);
  ud.expandDual(dt, kkt_matrix_ref, kkt_residual_ref, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(parnmpc.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(parnmpc.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  parnmpc.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


TEST_F(SplitUnconstrParNMPCTest, computeKKTResidual) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitUnconstrParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  parnmpc.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  parnmpc.computeKKTResidual(robot, t, dt, s_prev.q, s_prev.v, s, s_next, kkt_matrix, kkt_residual);
  SplitKKTResidual kkt_residual_ref(robot);
  SplitKKTMatrix kkt_matrix_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost = cost->linearizeStageCost(robot, cost_data, t, dt, s, kkt_residual_ref);
  constraints->linearizeConstraints(robot, constraints_data, dt, s, kkt_residual_ref);
  stage_cost += dt * constraints_data.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(dt, s_prev.q, s_prev.v, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  UnconstrDynamics ud(robot);
  ud.linearizeUnconstrDynamics(robot, dt, s, kkt_residual_ref);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual_ref.KKTError();
  kkt_error_ref += (dt*dt) * ud.KKTError();
  kkt_error_ref += (dt*dt) * constraints_data.KKTError();
  kkt_residual_ref.kkt_error = kkt_error_ref;
  EXPECT_DOUBLE_EQ(kkt_error_ref, parnmpc.KKTError(kkt_residual, dt));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitUnconstrParNMPCTest, evalOCP) {
  const auto s_prev = SplitSolution::Random(robot);
  const auto s = SplitSolution::Random(robot);
  const auto d = SplitDirection::Random(robot);
  SplitUnconstrParNMPC parnmpc(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  parnmpc.initConstraints(robot, 10, s);
  SplitKKTResidual kkt_residual(robot);
  parnmpc.evalOCP(robot, t, dt, s_prev.q, s_prev.v, s, kkt_residual);
  const double stage_cost = parnmpc.stageCost();
  const double constraint_violation = parnmpc.constraintViolation(kkt_residual, dt);
  SplitKKTResidual kkt_residual_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = cost->evalStageCost(robot, cost_data, t, dt, s);
  constraints->evalConstraint(robot, constraints_data, s);
  stage_cost_ref += dt * constraints_data.logBarrier();
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  unconstr::stateequation::computeBackwardEulerResidual(dt, s_prev.q, s_prev.v,
                                                        s, kkt_residual_ref);
  UnconstrDynamics ud(robot);
  ud.evalUnconstrDynamics(robot, s);
  double constraint_violation_ref = 0;
  constraint_violation_ref += kkt_residual_ref.constraintViolation();
  constraint_violation_ref += dt * constraints_data.constraintViolation();
  constraint_violation_ref += dt * ud.constraintViolation();
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}