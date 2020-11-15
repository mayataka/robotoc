#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/constraints/impulse_normal_force.hpp"


namespace idocp {

class SplitImpulseOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static std::shared_ptr<ImpulseCostFunction> createCost(const Robot& robot);

  static std::shared_ptr<ImpulseConstraints> createConstraints(const Robot& robot);

  static ImpulseSplitSolution generateFeasibleSolution(
      Robot& robot, const ImpulseStatus& impulse_status,
      const std::shared_ptr<ImpulseConstraints>& constraints);

  static void testLinearizeOCPAndRiccatiRecursion(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<ImpulseCostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints);

  static void testComputeKKTResidualEmptyCostAndEmptyConstraints(
      Robot& robot, const ImpulseStatus& impulse_status);

  static void testComputeKKTResidualEmptyCost(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<ImpulseConstraints>& constraints);

  static void testComputeKKTResidualEmptyConstraints(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<ImpulseCostFunction>& cost);

  static void testComputeKKTResidual(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<ImpulseCostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints);

  static void testCostAndConstraintViolation(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const std::shared_ptr<ImpulseCostFunction>& cost,
      const std::shared_ptr<ImpulseConstraints>& constraints);

  std::string fixed_base_urdf, floating_base_urdf;
};


std::shared_ptr<ImpulseCostFunction> SplitImpulseOCPTest::createCost(const Robot& robot) {
  auto joint_cost = std::make_shared<JointSpaceImpulseCost>(robot);
  const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_ref);
  const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
  const Eigen::VectorXd dv_ref = Eigen::VectorXd::Random(robot.dimv());
  joint_cost->set_q_weight(q_weight);
  joint_cost->set_q_ref(q_ref);
  joint_cost->set_v_weight(v_weight);
  joint_cost->set_v_ref(v_ref);
  joint_cost->set_dv_weight(dv_weight);
  joint_cost->set_dv_ref(dv_ref);
  auto impulse_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
  std::vector<Eigen::Vector3d> f_weight;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    f_weight.push_back(Eigen::Vector3d::Constant(0.001));
  }
  impulse_force_cost->set_f_weight(f_weight);
  auto cost = std::make_shared<ImpulseCostFunction>();
  cost->push_back(joint_cost);
  cost->push_back(impulse_force_cost);
  return cost;
}


std::shared_ptr<ImpulseConstraints> SplitImpulseOCPTest::createConstraints(const Robot& robot) {
  auto impulse_normal_force = std::make_shared<ImpulseNormalForce>(robot);
  auto constraints = std::make_shared<ImpulseConstraints>();
  constraints->push_back(impulse_normal_force);
  return constraints;
}


ImpulseSplitSolution SplitImpulseOCPTest::generateFeasibleSolution(
    Robot& robot, const ImpulseStatus& impulse_status,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  auto data = constraints->createConstraintsData(robot);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  while (!constraints->isFeasible(robot, data, s)) {
    s = ImpulseSplitSolution::Random(robot, impulse_status);
  }
  return s;
}


void SplitImpulseOCPTest::testLinearizeOCPAndRiccatiRecursion(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  RiccatiFactorization riccati_next(robot), riccati(robot), riccati_ref(robot);
  ImpulseSplitDirection d(robot), d_ref(robot);
  SplitDirection d_next(robot), d_next_ref(robot);
  d.setImpulseStatus(impulse_status);
  d_ref.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  const Eigen::MatrixXd P = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(2*dimv, 2*dimv);
  riccati_next.Pqq = P.topLeftCorner(dimv, dimv);
  riccati_next.Pqv = P.topRightCorner(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  riccati_next.Pvv = P.bottomRightCorner(dimv, dimv);
  ocp.linearizeOCP(robot, impulse_status, t, s_prev.q, s, s_next);
  ocp.backwardRiccatiRecursion(riccati_next, riccati);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  cost->computeStageCostHessian(robot, cost_data, t, s, kkt_matrix);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  constraints->condenseSlackAndDual(robot, constraints_data, s, kkt_matrix, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  id.condenseImpulseDynamics(robot, impulse_status, kkt_matrix, kkt_residual);
  ImpulseRiccatiFactorizer riccati_factorizer(robot);
  riccati_factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  d.dx().setRandom();
  riccati.Pi.setRandom();
  riccati.pi.setRandom();
  riccati.N.setRandom();
  riccati.n.setRandom();
  RiccatiFactorization riccati_next_ref = riccati_next;
  ocp.forwardRiccatiRecursionSerial(riccati, riccati_next);
  riccati_factorizer.forwardRiccatiRecursionSerial(riccati, kkt_matrix, kkt_residual, riccati_next_ref);
  EXPECT_TRUE(riccati_next.Pi.isApprox(riccati_next_ref.Pi));
  EXPECT_TRUE(riccati_next.pi.isApprox(riccati_next_ref.pi));
  EXPECT_TRUE(riccati_next.N.isApprox(riccati_next_ref.N));
  Eigen::MatrixXd Eq = Eigen::MatrixXd::Zero(impulse_status.dimp(), robot.dimv());
  Eigen::VectorXd e = Eigen::VectorXd::Zero(impulse_status.dimp());
  ocp.getStateConstraintFactorization(Eq, e);
  EXPECT_TRUE(Eq.isApprox(kkt_matrix.Pq()));
  EXPECT_TRUE(e.isApprox(kkt_residual.P()));
  const Eigen::MatrixXd T_next = Eigen::MatrixXd::Random(2*robot.dimv(), impulse_status.dimp());
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(2*robot.dimv(), impulse_status.dimp());
  ocp.backwardStateConstraintFactorization(T_next, T);
  if (!robot.has_floating_base()) {
    kkt_matrix.Fqq().setIdentity();
  }
  const Eigen::MatrixXd T_ref = kkt_matrix.Fxx().transpose() * T_next;
  EXPECT_TRUE(T.isApprox(T_ref));
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*robot.dimv());
  d_ref = d;
  ocp.computePrimalDirection(robot, riccati, s, dx0, d);
  riccati_factorizer.computeStateDirection(riccati, dx0, d_ref);
  riccati_factorizer.computeCostateDirection(riccati, d_ref);
  id.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  ocp.computeDualDirection(robot, d_next, d);
  id.computeCondensedDualDirection(robot, kkt_matrix, kkt_residual, 
                                   d_next_ref.dgmm(), d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyCostAndEmptyConstraints(
    Robot& robot, const ImpulseStatus& impulse_status) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, empty_constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq());
    robot.dSubtractdConfigurationMinus(s_prev.q, s.q, kkt_matrix.Fqq_prev);
    kkt_residual.Fv() = s.v + s.dv - s_next.v;
    kkt_residual.lq() = kkt_matrix.Fqq().transpose() * s_next.lmd + kkt_matrix.Fqq_prev.transpose() * s.lmd;
    kkt_residual.lv() = s_next.gmm - s.gmm;
    kkt_residual.ldv = s_next.gmm;
  }
  else {
    kkt_residual.Fq() = s.q - s_next.q;
    kkt_residual.Fv() = s.v + s.dv - s_next.v;
    kkt_residual.lq() = s_next.lmd - s.lmd;
    kkt_residual.lv() = s_next.gmm - s.gmm;
    kkt_residual.ldv = s_next.gmm;
  }
  Eigen::VectorXd ImD = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dImD_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dImD_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, ImD);
  robot.RNEAImpulseDerivatives(s.q, s.dv, dImD_dq, dImD_ddv);
  kkt_residual.lq() += dImD_dq.transpose() * s.beta;
  kkt_residual.ldv += dImD_ddv.transpose() * s.beta;
  Eigen::VectorXd P = Eigen::VectorXd::Zero(impulse_status.dimp());
  Eigen::VectorXd V = Eigen::VectorXd::Zero(impulse_status.dimp());
  if (impulse_status.hasActiveImpulse()) {
    Eigen::MatrixXd dImD_df = Eigen::MatrixXd::Zero(robot.dimv(), impulse_status.dimp());
    robot.updateKinematics(s.q, s.v);
    ContactStatus contact_status(impulse_status.max_point_contacts());
    for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
      if (impulse_status.isImpulseActive(i)) {
        contact_status.activateContact(i);
      }
    }
    robot.dRNEAPartialdFext(contact_status, dImD_df);
    kkt_residual.lf() += dImD_df.transpose() * s.beta;
    Eigen::MatrixXd dV_dq = Eigen::MatrixXd::Zero(impulse_status.dimp(), robot.dimv());
    Eigen::MatrixXd dV_dv = Eigen::MatrixXd::Zero(impulse_status.dimp(), robot.dimv());
    Eigen::MatrixXd dP_dq = Eigen::MatrixXd::Zero(impulse_status.dimp(), robot.dimv());
    robot.updateKinematics(s.q, s.v);
    robot.computeImpulseVelocityResidual(impulse_status, V);
    robot.computeImpulseVelocityDerivatives(impulse_status, dV_dq, dV_dv);
    robot.computeImpulseConditionResidual(impulse_status, P);
    robot.computeImpulseConditionDerivative(impulse_status, dP_dq);
    kkt_residual.lq() += dV_dq.transpose() * s.mu_stack() + dP_dq.transpose() * s.xi_stack();
    kkt_residual.lv() += dV_dv.transpose() * s.mu_stack();
    kkt_residual.ldv += dV_dv.transpose() * s.mu_stack();
  }
  const double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                                + kkt_residual.lx().squaredNorm()
                                + kkt_residual.ldv.squaredNorm()
                                + kkt_residual.lf().squaredNorm()
                                + ImD.squaredNorm()
                                + V.squaredNorm()
                                + P.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyCost(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitImpulseOCP ocp(robot, empty_cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  const double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                                + kkt_residual.lx().squaredNorm()
                                + kkt_residual.ldv.squaredNorm()
                                + kkt_residual.lf().squaredNorm()
                                + id.squaredNormImpulseDynamicsResidual(kkt_residual)
                                + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyConstraints(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost) {
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, empty_constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitImpulseOCP ocp(robot, cost, empty_constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  robot.updateKinematics(s.q, s.v);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  const double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                                + kkt_residual.lx().squaredNorm()
                                + kkt_residual.ldv.squaredNorm()
                                + kkt_residual.lf().squaredNorm()
                                + id.squaredNormImpulseDynamicsResidual(kkt_residual);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testComputeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const SplitSolution s_next = SplitSolution::Random(robot);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  ocp.computeKKTResidual(robot, impulse_status, t, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamicsForwardEuler id(robot);
  robot.updateKinematics(s.q, s.v);
  id.linearizeImpulseDynamics(robot, impulse_status, s, kkt_matrix, kkt_residual);
  const double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                                + kkt_residual.lx().squaredNorm()
                                + kkt_residual.ldv.squaredNorm()
                                + kkt_residual.lf().squaredNorm()
                                + id.squaredNormImpulseDynamicsResidual(kkt_residual)
                                + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testCostAndConstraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  const SplitSolution s_prev = SplitSolution::Random(robot);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size = 0.3;
  ocp.initConstraints(robot, s);
  const double stage_cost = ocp.stageCost(robot, t, s, step_size);
  const double constraint_violation = ocp.constraintViolation(robot, impulse_status, t, s, s_prev.q, s_prev.v);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot);
  constraints->setSlackAndDual(robot, constraints_data, s);
  robot.updateKinematics(s.q, s.v);
  double stage_cost_ref = 0;
  stage_cost_ref += cost->l(robot, cost_data, t, s);
  stage_cost_ref += constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(stage_cost, stage_cost_ref);
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  stateequation::ComputeImpulseForwardEulerResidual(robot, s, s_prev.q, 
                                                    s_prev.v, kkt_residual);
  ImpulseDynamicsForwardEuler id(robot);
  id.computeImpulseDynamicsResidual(robot, impulse_status, s, kkt_residual);
  double constraint_violation_ref = 0;
  constraint_violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  constraint_violation_ref += stateequation::L1NormStateEuqationResidual(kkt_residual);
  constraint_violation_ref += id.l1NormImpulseDynamicsResidual(kkt_residual);
  EXPECT_DOUBLE_EQ(constraint_violation, constraint_violation_ref);
}


TEST_F(SplitImpulseOCPTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ImpulseStatus impulse_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCPAndRiccatiRecursion(robot, impulse_status, cost, constraints);
  testComputeKKTResidualEmptyCostAndEmptyConstraints(robot, impulse_status);
  testComputeKKTResidualEmptyCost(robot, impulse_status, constraints);
  testComputeKKTResidualEmptyConstraints(robot, impulse_status, cost);
  testComputeKKTResidual(robot, impulse_status, cost, constraints);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints);
  impulse_status.setImpulseStatus({true});
  testLinearizeOCPAndRiccatiRecursion(robot, impulse_status, cost, constraints);
  testComputeKKTResidualEmptyCostAndEmptyConstraints(robot, impulse_status);
  testComputeKKTResidualEmptyCost(robot, impulse_status, constraints);
  testComputeKKTResidualEmptyConstraints(robot, impulse_status, cost);
  testComputeKKTResidual(robot, impulse_status, cost, constraints);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints);
}


TEST_F(SplitImpulseOCPTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ImpulseStatus impulse_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false, false, false, false});
  const auto cost = createCost(robot);
  const auto constraints = createConstraints(robot);
  testLinearizeOCPAndRiccatiRecursion(robot, impulse_status, cost, constraints);
  testComputeKKTResidualEmptyCostAndEmptyConstraints(robot, impulse_status);
  testComputeKKTResidualEmptyCost(robot, impulse_status, constraints);
  testComputeKKTResidualEmptyConstraints(robot, impulse_status, cost);
  testComputeKKTResidual(robot, impulse_status, cost, constraints);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  impulse_status.setImpulseStatus(is_contact_active);
  testLinearizeOCPAndRiccatiRecursion(robot, impulse_status, cost, constraints);
  testComputeKKTResidualEmptyCostAndEmptyConstraints(robot, impulse_status);
  testComputeKKTResidualEmptyCost(robot, impulse_status, constraints);
  testComputeKKTResidualEmptyConstraints(robot, impulse_status, cost);
  testComputeKKTResidual(robot, impulse_status, cost, constraints);
  testCostAndConstraintViolation(robot, impulse_status, cost, constraints);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}