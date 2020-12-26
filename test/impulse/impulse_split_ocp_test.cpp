#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/constraints/impulse_normal_force.hpp"


namespace idocp {

class ImpulseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);

  static std::shared_ptr<ImpulseConstraints> createConstraints(const Robot& robot);

  static ImpulseSplitSolution generateFeasibleSolution(
      Robot& robot, const ImpulseStatus& impulse_status,
      const std::shared_ptr<ImpulseConstraints>& constraints);

  static void testLinearizeOCP(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints, 
      const bool is_state_constraint_valid);

  static void testComputeKKTResidual(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints,
      const bool is_state_constraint_valid);

  static void testCostAndConstraintViolation(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints,
      const bool is_state_constraint_valid);

  std::string fixed_base_urdf, floating_base_urdf;
};


std::shared_ptr<CostFunction> ImpulseSplitOCPTest::createCost(const Robot& robot) {
  auto config_cost = std::make_shared<ConfigurationSpaceCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_ref(v_ref);
  config_cost->set_dvi_weight(dv_weight);
  auto contact_force_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_force_cost->set_f_weight(f_weight);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(config_cost);
  cost->push_back(contact_force_cost);
  return cost;
}


std::shared_ptr<ImpulseConstraints> ImpulseSplitOCPTest::createConstraints(const Robot& robot) {
  auto impulse_normal_force = std::make_shared<ImpulseNormalForce>(robot);
  auto constraints = std::make_shared<ImpulseConstraints>();
  constraints->push_back(impulse_normal_force);
  return constraints;
}


ImpulseSplitSolution ImpulseSplitOCPTest::generateFeasibleSolution(
    Robot& robot, const ImpulseStatus& impulse_status,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  auto data = constraints->createConstraintsData(robot);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  while (!constraints->isFeasible(robot, data, s)) {
    s = ImpulseSplitSolution::Random(robot, impulse_status);
  }
  return s;
}


void ImpulseSplitOCPTest::testLinearizeOCP(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints,
    const bool is_state_constraint_valid) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ocp.linearizeOCP(robot, impulse_status, t, s_prev.q, s, s_next, kkt_matrix, kkt_residual, is_state_constraint_valid);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  cost->computeImpulseCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  cost->computeImpulseCostHessian(robot, cost_data, t, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref, is_state_constraint_valid );
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  ocp.computeCondensedPrimalDirection(robot, s, d);
  id.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, d_next, d);
  id.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, d_next.dgmm(), d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated, is_state_constraint_valid);
  s_updated_ref.integrate(robot, step_size, d, is_state_constraint_valid);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void ImpulseSplitOCPTest::testComputeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints, 
    const bool is_state_constraint_valid) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next, kkt_matrix, kkt_residual, is_state_constraint_valid);
  const double kkt_error = ocp.squaredNormKKTResidual(kkt_residual, is_state_constraint_valid);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  cost->computeImpulseCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual_ref);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, v_after_impulse);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix_ref, kkt_residual_ref, is_state_constraint_valid);
  const double kkt_error_ref = kkt_residual_ref.Fx().squaredNorm()
                                + kkt_residual_ref.lx().squaredNorm()
                                + kkt_residual_ref.ldv.squaredNorm()
                                + kkt_residual_ref.lf().squaredNorm()
                                + id.squaredNormImpulseDynamicsResidual(kkt_residual_ref, is_state_constraint_valid)
                                + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void ImpulseSplitOCPTest::testCostAndConstraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints,
    const bool is_state_constraint_valid) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  ImpulseSplitOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size = 0.3;
  ocp.initConstraints(robot, s);
  const double stage_cost = ocp.stageCost(robot, t, s, step_size);
  ImpulseSplitKKTResidual kkt_residual(robot);
  const double constraint_violation = ocp.constraintViolation(robot, impulse_status, t, s, s_prev.q, s_prev.v, kkt_residual, is_state_constraint_valid);
  ImpulseSplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setImpulseStatus(impulse_status);
  ImpulseSplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const Eigen::VectorXd v_after_impulse = s.v + s.dv;
  robot.updateKinematics(s.q, v_after_impulse);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeImpulseCost(robot, cost_data, t, s);
  stage_cost_ref += constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::ComputeImpulseForwardEulerResidual(robot, s, s_prev.q, 
                                                    s_prev.v, kkt_residual_ref);
  ImpulseDynamicsForwardEuler id(robot);
  id.computeImpulseDynamicsResidual(robot, impulse_status, s, kkt_residual_ref, is_state_constraint_valid);
  double constraint_violation_ref = 0;
  constraint_violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::L1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += id.l1NormImpulseDynamicsResidual(kkt_residual_ref, is_state_constraint_valid);
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}


TEST_F(ImpulseSplitOCPTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ImpulseStatus impulse_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    impulse_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(fixed_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, impulse_status, cost, constraints, false);
  testLinearizeOCP(robot, impulse_status, cost, constraints, true);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, false);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, false);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, true);
  impulse_status.setImpulseStatus({true});
  testLinearizeOCP(robot, impulse_status, cost, constraints, false);
  testLinearizeOCP(robot, impulse_status, cost, constraints, true);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, false);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, false);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, true);
}


TEST_F(ImpulseSplitOCPTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ImpulseStatus impulse_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    impulse_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(floating_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false, false, false, false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, impulse_status, cost, constraints, false);
  testLinearizeOCP(robot, impulse_status, cost, constraints, true);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, false);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, false);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, true);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  impulse_status.setImpulseStatus(is_contact_active);
  testLinearizeOCP(robot, impulse_status, cost, constraints, false);
  testLinearizeOCP(robot, impulse_status, cost, constraints, true);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, false);
  testComputeKKTResidual(robot, impulse_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, false);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints, true);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}