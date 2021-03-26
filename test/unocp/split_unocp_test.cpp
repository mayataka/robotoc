#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unocp.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class SplitUnOCPTest : public ::testing::Test {
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


TEST_F(SplitUnOCPTest, linearizeOCP) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitUnOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitUnKKTMatrix unkkt_matrix(robot);
  SplitUnKKTResidual unkkt_residual(robot);
  ocp.linearizeOCP(robot, t, dt, s_prev.q, s, s_next, unkkt_matrix, unkkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  SplitUnKKTMatrix unkkt_matrix_ref(robot);
  SplitUnKKTResidual unkkt_residual_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  cost->computeStageCostHessian(robot, cost_data, t, dt, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dt, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  UnconstrainedDynamics ud(robot);
  ud.linearizeUnconstrainedDynamics(robot, dt, s, kkt_residual_ref);
  ud.condenseUnconstrainedDynamics(kkt_matrix_ref, kkt_residual_ref, unkkt_matrix_ref, unkkt_residual_ref);
  EXPECT_TRUE(unkkt_matrix.isApprox(unkkt_matrix_ref));
  EXPECT_TRUE(unkkt_residual.isApprox(unkkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot);
  auto d_ref = d;
  ocp.computeCondensedDirection(robot, dt, s, d);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  ud.computeCondensedDirection(dt, kkt_matrix_ref, kkt_residual_ref, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


TEST_F(SplitUnOCPTest, computeKKTResidual) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitUnOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  ocp.computeKKTResidual(robot, t, dt, s_prev.q, s, s_next);
  SplitKKTResidual kkt_residual_ref(robot);
  SplitKKTMatrix kkt_matrix_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dt, s, kkt_residual_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, dt, s, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dt, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  UnconstrainedDynamics ud(robot);
  ud.linearizeUnconstrainedDynamics(robot, dt, s, kkt_residual_ref);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual_ref.lx().squaredNorm();
  kkt_error_ref += kkt_residual_ref.la.squaredNorm();
  kkt_error_ref += kkt_residual_ref.lu().squaredNorm();
  kkt_error_ref += stateequation::squaredNormStateEuqationResidual(kkt_residual_ref);
  kkt_error_ref += ud.squaredNormUnconstrainedDynamicsResidual(dt);
  kkt_error_ref += dt * dt * constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error_ref, ocp.squaredNormKKTResidual(dt));
}


TEST_F(SplitUnOCPTest, costAndConstraintViolation) {
  const SplitSolution s = SplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  const SplitDirection d = SplitDirection::Random(robot);
  SplitUnOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size = 0.3;
  ocp.initConstraints(robot, 10, s);
  const double stage_cost = ocp.stageCost(robot, t, dt, s, step_size);
  const double constraint_violation = ocp.constraintViolation(robot, t, dt, s, s_next.q, s_next.v);
  SplitKKTResidual kkt_residual_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeStageCost(robot, cost_data, t, dt, s);
  stage_cost_ref += dt * constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::computeForwardEulerResidual(robot, dt, s, s_next.q, 
                                             s_next.v, kkt_residual_ref);
  UnconstrainedDynamics cd(robot);
  cd.computeUnconstrainedDynamicsResidual(robot, s);
  double constraint_violation_ref = 0;
  constraint_violation_ref += dt * constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += cd.l1NormUnconstrainedDynamicsResidual(dt);
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}