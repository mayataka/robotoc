#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/terminal_parnmpc.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/configuration_space_cost.hpp"
#include "idocp/cost/task_space_3d_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"


namespace idocp {

class TerminalParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    baumgarte_time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);

  static std::shared_ptr<Constraints> createConstraints(const Robot& robot);

  static SplitSolution generateFeasibleSolution(
      Robot& robot, const ContactStatus& contact_staus,
      const std::shared_ptr<Constraints>& constraints);

  void testLinearizeOCP(
      Robot& robot, const ContactStatus& contact_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints) const;

  void testComputeKKTResidual(
      Robot& robot, const ContactStatus& contact_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints) const;

  void testCostAndConstraintViolation(
      Robot& robot, const ContactStatus& contact_status, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double baumgarte_time_step;
};


std::shared_ptr<CostFunction> TerminalParNMPCTest::createCost(const Robot& robot) {
  auto config_cost = std::make_shared<ConfigurationSpaceCost >(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimu()).array().abs();
  const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  config_cost->set_q_weight(q_weight);
  config_cost->set_q_ref(q_ref);
  config_cost->set_v_weight(v_weight);
  config_cost->set_v_ref(v_ref);
  config_cost->set_a_weight(a_weight);
  config_cost->set_u_weight(u_weight);
  config_cost->set_u_ref(u_ref);
  config_cost->set_qf_weight(qf_weight);
  config_cost->set_vf_weight(vf_weight);
  const int task_frame = 10;
  auto task_space_3d_cost = std::make_shared<TaskSpace3DCost >(robot, task_frame);
  const Eigen::Vector3d q_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d qf_3d_weight = Eigen::Vector3d::Random().array().abs();
  const Eigen::Vector3d q_3d_ref = Eigen::Vector3d::Random();
  task_space_3d_cost->set_q_3d_weight(q_3d_weight);
  task_space_3d_cost->set_qf_3d_weight(qf_3d_weight);
  task_space_3d_cost->set_q_3d_ref(q_3d_ref);
  auto contact_force_cost = std::make_shared<idocp::ContactForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<robot.maxPointContacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  contact_force_cost->set_f_weight(f_weight);
  auto cost = std::make_shared<CostFunction>();
  cost->push_back(config_cost);
  cost->push_back(task_space_3d_cost);
  cost->push_back(contact_force_cost);
  return cost;
}


std::shared_ptr<Constraints> TerminalParNMPCTest::createConstraints(const Robot& robot) {
  auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
  auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
  auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
  auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
  auto constraints = std::make_shared<Constraints>();
  constraints->push_back(joint_upper_limit); 
  constraints->push_back(joint_lower_limit);
  constraints->push_back(velocity_lower_limit); 
  constraints->push_back(velocity_upper_limit);
  return constraints;
}


SplitSolution TerminalParNMPCTest::generateFeasibleSolution(
    Robot& robot, const ContactStatus& contact_staus,
    const std::shared_ptr<Constraints>& constraints) {
  auto data = constraints->createConstraintsData(robot, 10);
  SplitSolution s = SplitSolution::Random(robot, contact_staus);
  while (!constraints->isFeasible(robot, data, s)) {
    s = SplitSolution::Random(robot, contact_staus);
  }
  return s;
}


void TerminalParNMPCTest::testLinearizeOCP(
    Robot& robot, const ContactStatus& contact_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints) const {
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = generateFeasibleSolution(robot, contact_status, constraints);
  TerminalParNMPC ocp(robot, cost, constraints, baumgarte_time_step);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.linearizeOCP(robot, contact_status, t, dtau, s_prev.q, s_prev.v, s, kkt_matrix, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual_ref);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix_ref);
  cost->computeTerminalCostHessian(robot, cost_data, t, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeBackwardEulerTerminal(robot, dtau, s_prev.q, s_prev.v, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseBackwardEuler(robot, dtau, s_prev.q, s, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual_ref);
  cd.condenseContactDynamicsBackwardEuler(robot, contact_status, dtau, kkt_matrix_ref, kkt_residual_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  SplitDirection d = SplitDirection::Random(robot, contact_status);
  auto d_ref = d;
  ocp.computeCondensedPrimalDirection(robot, dtau, s, d);
  cd.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, d);
  cd.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, d.dgmm(), d_ref);
  Eigen::VectorXd dlmd_ref = d_ref.dlmd();
  if (robot.hasFloatingBase()) {
    d_ref.dlmd().head(6) = kkt_matrix_ref.Fqq_prev_inv.transpose() * dlmd_ref.head(6);
  }
  EXPECT_TRUE(d.isApprox(d_ref));
  const double step_size = std::abs(Eigen::VectorXd::Random(1)[0]);
  auto s_updated = s;
  auto s_updated_ref = s;
  ocp.updatePrimal(robot, step_size, d, s_updated);
  s_updated_ref.integrate(robot, step_size, d);
  constraints->updateSlack(constraints_data, step_size);
  EXPECT_TRUE(s_updated.isApprox(s_updated_ref));
}


void TerminalParNMPCTest::testComputeKKTResidual(
    Robot& robot, const ContactStatus& contact_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints) const {
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  TerminalParNMPC ocp(robot, cost, constraints, baumgarte_time_step);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, s_prev.q, s_prev.v, s, kkt_matrix, kkt_residual);
  const double kkt_error = ocp.squaredNormKKTResidual(kkt_residual, dtau);
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual_ref);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual_ref);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual_ref);
  stateequation::linearizeBackwardEulerTerminal(robot, dtau, s_prev.q, s_prev.v, s, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual_ref);
  double kkt_error_ref = kkt_residual_ref.Fx().squaredNorm()
                         + kkt_residual_ref.lx().squaredNorm()
                         + kkt_residual_ref.la.squaredNorm()
                         + kkt_residual_ref.lf().squaredNorm()
                         + kkt_residual_ref.lu().squaredNorm()
                         + cd.squaredNormContactDynamicsResidual(dtau)
                         + dtau * dtau * constraints->squaredNormPrimalAndDualResidual(constraints_data);
  if (robot.hasFloatingBase()) {
    kkt_error_ref += kkt_residual_ref.lu_passive.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void TerminalParNMPCTest::testCostAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints) const {
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
  TerminalParNMPC ocp(robot, cost, constraints, baumgarte_time_step);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size = 0.3;
  ocp.initConstraints(robot, 10, s);
  const double stage_cost = ocp.stageCost(robot, t, dtau, s, step_size);
  SplitKKTResidual kkt_residual(robot);
  const double constraint_violation = ocp.constraintViolation(robot, contact_status, t, dtau, s_prev.q, s_prev.v, s, kkt_residual);
  SplitKKTMatrix kkt_matrix_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_residual_ref.setContactStatus(contact_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->computeStageCost(robot, cost_data, t, dtau, s);
  stage_cost_ref += dtau * constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::computeBackwardEulerResidual(robot, dtau, s_prev.q, s_prev.v, s, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  cd.computeContactDynamicsResidual(robot, contact_status, s);
  double constraint_violation_ref = 0;
  constraint_violation_ref += dtau * constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += cd.l1NormContactDynamicsResidual(dtau);
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(TerminalParNMPCTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(fixed_base_urdf, contact_frames);
  contact_status.setContactStatus({false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  contact_status.setContactStatus({true});
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
}


TEST_F(TerminalParNMPCTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  for (int i=0; i<contact_frames.size(); ++i) {
    contact_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  Robot robot(floating_base_urdf, contact_frames);
  contact_status.setContactStatus({false, false, false, false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  contact_status.setContactStatus(is_contact_active);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}