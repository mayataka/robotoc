#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_ocp.hpp"
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
#include "idocp/ocp/state_equation.hpp"
#include "idocp/ocp/contact_dynamics.hpp"
#include "idocp/ocp/forward_switching_constraint.hpp"


namespace idocp {

class SplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    baumgarte_time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau_next = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<CostFunction> createCost(const Robot& robot);

  static std::shared_ptr<Constraints> createConstraints(const Robot& robot);

  static SplitSolution generateFeasibleSolution(
      Robot& robot, const ContactStatus& contact_staus,
      const std::shared_ptr<Constraints>& constraints);

  static SplitSolution generateFeasibleSolution(
      Robot& robot, const ContactStatus& contact_staus, 
      const ImpulseStatus& impulse_status,
      const std::shared_ptr<Constraints>& constraints);

  void testLinearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                        const std::shared_ptr<CostFunction>& cost,
                        const std::shared_ptr<Constraints>& constraints,
                        const bool switching_constraint=false) const;

  void testComputeKKTResidual(Robot& robot, const ContactStatus& contact_status, 
                              const std::shared_ptr<CostFunction>& cost,
                              const std::shared_ptr<Constraints>& constraints,
                              const bool switching_constraint=false) const;

  void testCostAndConstraintViolation(Robot& robot, const ContactStatus& contact_status, 
                                      const std::shared_ptr<CostFunction>& cost,
                                      const std::shared_ptr<Constraints>& constraints,
                                      const bool switching_constraint=false) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double baumgarte_time_step, t, dtau, dtau_next;
};


std::shared_ptr<CostFunction> SplitOCPTest::createCost(const Robot& robot) {
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


std::shared_ptr<Constraints> SplitOCPTest::createConstraints(const Robot& robot) {
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


SplitSolution SplitOCPTest::generateFeasibleSolution(
    Robot& robot, const ContactStatus& contact_staus,
    const std::shared_ptr<Constraints>& constraints) {
  auto data = constraints->createConstraintsData(robot, 10);
  SplitSolution s = SplitSolution::Random(robot, contact_staus);
  while (!constraints->isFeasible(robot, data, s)) {
    s = SplitSolution::Random(robot, contact_staus);
  }
  return s;
}


SplitSolution SplitOCPTest::generateFeasibleSolution(
    Robot& robot, const ContactStatus& contact_staus,
    const ImpulseStatus& impulse_status,
    const std::shared_ptr<Constraints>& constraints) {
  auto data = constraints->createConstraintsData(robot, 10);
  SplitSolution s = SplitSolution::Random(robot, contact_staus, impulse_status);
  while (!constraints->isFeasible(robot, data, s)) {
    s = SplitSolution::Random(robot, contact_staus);
  }
  return s;
}


void SplitOCPTest::testLinearizeOCP(Robot& robot, const ContactStatus& contact_status, 
                                    const std::shared_ptr<CostFunction>& cost,
                                    const std::shared_ptr<Constraints>& constraints,
                                    const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = generateFeasibleSolution(robot, contact_status, impulse_status, constraints);
  }
  else {
    stmp = generateFeasibleSolution(robot, contact_status, constraints);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  SplitOCP ocp(robot, cost, constraints, baumgarte_time_step);
  ocp.initConstraints(robot, 10, s);
  const int dimv = robot.dimv();
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitStateConstraintJacobian jac(robot);
  if (switching_constraint) {
    ocp.linearizeOCP(robot, contact_status, t, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual, 
                     impulse_status, dtau_next, jac);
  }
  else {
    ocp.linearizeOCP(robot, contact_status, t, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  }
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  SplitStateConstraintJacobian jac_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual_ref);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix_ref);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual_ref);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, kkt_matrix_ref, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  stateequation::condenseForwardEuler(robot, dtau, s, s_next.q, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual_ref);
  if (switching_constraint) {
    ForwardSwitchingConstraint sc(robot);
    sc.linearizeSwitchingConstraint(robot, impulse_status, dtau, dtau_next, s, 
                                    kkt_matrix_ref, kkt_residual_ref, jac_ref);
    cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix_ref, kkt_residual_ref, true);
    cd.condenseSwitchingConstraint(kkt_residual_ref, jac_ref);
    EXPECT_TRUE(jac.isApprox(jac_ref));
  }
  else {
    cd.condenseContactDynamics(robot, contact_status, dtau, kkt_matrix_ref, kkt_residual_ref, true);
  }
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_residual.lx().isApprox(kkt_residual_ref.lx()));
  EXPECT_TRUE(kkt_residual.lu().isApprox(kkt_residual_ref.lu()));
  EXPECT_TRUE(kkt_residual.la.isApprox(kkt_residual_ref.la));
  EXPECT_TRUE(kkt_residual.P().isApprox(kkt_residual_ref.P()));
  SplitDirection d;
  if (switching_constraint) {
    d = SplitDirection::Random(robot, contact_status, impulse_status);
  }
  else {
    d = SplitDirection::Random(robot, contact_status);
  }
  auto d_ref = d;
  const SplitDirection d_next = SplitDirection::Random(robot);
  ocp.computeCondensedPrimalDirection(robot, dtau, s, d);
  cd.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  EXPECT_DOUBLE_EQ(ocp.maxPrimalStepSize(), constraints->maxSlackStepSize(constraints_data));
  EXPECT_DOUBLE_EQ(ocp.maxDualStepSize(), constraints->maxDualStepSize(constraints_data));
  ocp.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, d_next, d);
  cd.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, d_next.dgmm(), d_ref);
  Eigen::VectorXd dlmd_ref = d_ref.dlmd();
  if (robot.hasFloatingBase()) {
    d_ref.dlmd().head(6) = - kkt_matrix_ref.Fqq_prev_inv.transpose() * dlmd_ref.head(6);
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


void SplitOCPTest::testComputeKKTResidual(
    Robot& robot, const ContactStatus& contact_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints,
    const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = generateFeasibleSolution(robot, contact_status, impulse_status, constraints);
  }
  else {
    stmp = generateFeasibleSolution(robot, contact_status, constraints);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  SplitOCP ocp(robot, cost, constraints, baumgarte_time_step);
  ocp.initConstraints(robot, 10, s);
  SplitKKTMatrix kkt_matrix(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitStateConstraintJacobian jac(robot);
  if (switching_constraint) {
    ocp.computeKKTResidual(robot, contact_status, t, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual, 
                           impulse_status, dtau_next, jac);
  }
  else {
    ocp.computeKKTResidual(robot, contact_status, t, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  }
  const double kkt_error = ocp.squaredNormKKTResidual(kkt_residual, dtau);
  SplitKKTMatrix kkt_matrix_ref(robot);
  SplitKKTResidual kkt_residual_ref(robot);
  kkt_matrix_ref.setContactStatus(contact_status);
  kkt_residual_ref.setContactStatus(contact_status);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
  }
  SplitStateConstraintJacobian jac_ref(robot);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual_ref);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual_ref);
  stateequation::linearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix_ref, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_residual_ref);
  if (switching_constraint) {
    kkt_matrix_ref.setImpulseStatus(impulse_status);
    kkt_residual_ref.setImpulseStatus(impulse_status);
    ForwardSwitchingConstraint sc(robot);
    sc.linearizeSwitchingConstraint(robot, impulse_status, dtau, dtau_next, s, 
                                    kkt_matrix_ref, kkt_residual_ref, jac_ref);
    EXPECT_TRUE(jac.isApprox(jac_ref));
  }
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
  if (switching_constraint) {
    kkt_error_ref += kkt_residual_ref.P().squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitOCPTest::testCostAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints,
    const bool switching_constraint) const {
  ImpulseStatus impulse_status;
  if (switching_constraint) {
    impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    if (!impulse_status.hasActiveImpulse()) {
      impulse_status.activateImpulse(0);
    }
  }
  SplitSolution stmp;
  if (switching_constraint) {
    stmp = generateFeasibleSolution(robot, contact_status, impulse_status, constraints);
  }
  else {
    stmp = generateFeasibleSolution(robot, contact_status, constraints);
  }
  const SplitSolution s_prev = SplitSolution::Random(robot, contact_status);
  const SplitSolution s = stmp;
  const SplitSolution s_next = SplitSolution::Random(robot, contact_status);
  SplitOCP ocp(robot, cost, constraints, baumgarte_time_step);
  const double step_size = 0.3;
  ocp.initConstraints(robot, 10, s);
  const double stage_cost = ocp.stageCost(robot, t, dtau, s, step_size);
  SplitKKTResidual kkt_residual(robot);
  double constraint_violation;
  if (switching_constraint) {
    constraint_violation = ocp.constraintViolation(robot, contact_status, t, dtau, s, s_next.q, s_next.v, kkt_residual, impulse_status, dtau_next);
  }
  else {
    constraint_violation = ocp.constraintViolation(robot, contact_status, t, dtau, s, s_next.q, s_next.v, kkt_residual);
  }
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
  stateequation::computeForwardEulerResidual(robot, dtau, s, s_next.q, 
                                             s_next.v, kkt_residual_ref);
  ContactDynamics cd(robot, baumgarte_time_step);
  cd.computeContactDynamicsResidual(robot, contact_status, s);
  if (switching_constraint) {
    kkt_residual_ref.setImpulseStatus(impulse_status);
    ForwardSwitchingConstraint sc(robot);
    sc.computeSwitchingConstraintResidual(robot, impulse_status, dtau, dtau_next, 
                                          s, kkt_residual_ref);
  }
  double constraint_violation_ref = 0;
  constraint_violation_ref += dtau * constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::l1NormStateEuqationResidual(kkt_residual_ref);
  constraint_violation_ref += cd.l1NormContactDynamicsResidual(dtau);
  if (switching_constraint) {
    constraint_violation_ref += kkt_residual_ref.P().lpNorm<1>();
  }
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitOCPTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  auto contact_status = robot.createContactStatus();
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  contact_status.activateContact(0);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  testLinearizeOCP(robot, contact_status, cost, constraints, true);
  testComputeKKTResidual(robot, contact_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints, true);
}


TEST_F(SplitOCPTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  ContactStatus contact_status = robot.createContactStatus();
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  testLinearizeOCP(robot, contact_status, cost, constraints);
  testComputeKKTResidual(robot, contact_status, cost, constraints);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints);
  testLinearizeOCP(robot, contact_status, cost, constraints, true);
  testComputeKKTResidual(robot, contact_status, cost, constraints, true);
  testCostAndConstraintViolation(robot, contact_status, cost, constraints, true);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}