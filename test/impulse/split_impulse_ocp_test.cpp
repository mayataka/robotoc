#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/ocp/riccati_solution.hpp"
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

  static SplitSolution generateFeasibleSolution(
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
  joint_cost->set_dv_weight(a_weight);
  joint_cost->set_dv_ref(a_ref);
  auto impusle_force_cost = std::make_shared<idocp::ImpulseForceCost>(robot);
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


SplitSolution SplitImpulseOCPTest::generateFeasibleSolution(
    Robot& robot, const ImpulseStatus& impulse_status,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  SplitImpulseOCP ocp(robot, empty_cost, constraints);
  ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  while (!ocp.isFeasible(robot, s)) {
    s = ImpulseSplitSolution::Random(robot, impulse_status);
  }
  return s;
}


void SplitImpulseOCPTest::testLinearizeOCPAndRiccatiRecursion(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<CostFunction>& cost,
    const std::shared_ptr<Constraints>& constraints) {
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = generateFeasibleSolution(robot, impulse_status, constraints);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, s);
  RiccatiSolution riccati_next(robot), riccati(robot), riccati_ref(robot);
  ImpulseSplitDirection
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
  ocp.backwardRiccatiRecursion(dtau, riccati_next, riccati);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, kkt_matrix, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, impulse_status, dtau, s, kkt_matrix, kkt_residual);
  cd.condenseContactDynamics(robot, impulse_status, dtau, kkt_matrix, kkt_residual);
  RiccatiFactorizer riccati_factorizer(robot);
  RiccatiGain riccati_gain(robot);
  riccati_factorizer.backwardRiccatiRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, riccati_gain, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  d.dx().setRandom();
  d_ref = d;
  ocp.forwardRiccatiRecursion(dtau, d, d_next);
  riccati_factorizer.computeControlInputDirection(riccati_gain, d_ref);
  riccati_factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d_ref, dtau, d_next_ref);
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
  ocp.computeCondensedPrimalDirection(robot, dtau, riccati_ref, s, d);
  riccati_factorizer.computeCostateDirection(riccati, d_ref);
  cd.computeCondensedPrimalDirection(robot, d_ref);
  constraints->computeSlackAndDualDirection(robot, constraints_data, dtau, s, d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
  ocp.computeCondensedDualDirection(robot, dtau, d_next, d);
  cd.computeCondensedDualDirection(robot, dtau, kkt_matrix, kkt_residual, 
                                   d_next_ref.dgmm(), d_ref);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyCostAndEmptyConstraints(
    Robot& robot, const ImpulseStatus& impulse_status) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, dtau, s);
  ocp.computeKKTResidual(robot, impulse_status, t, dtau, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
    robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq());
    robot.dSubtractdConfigurationMinus(s_prev.q, s.q, kkt_matrix.Fqq_prev);
    kkt_residual.Fq() += dtau * s.v;
    kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
    kkt_residual.lq() = kkt_matrix.Fqq().transpose() * s_next.lmd + kkt_matrix.Fqq_prev.transpose() * s.lmd;
    kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
    kkt_residual.la = dtau * s_next.gmm;
  }
  else {
    kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
    kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
    kkt_residual.lq() = s_next.lmd - s.lmd;
    kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
    kkt_residual.la = dtau * s_next.gmm;
  }
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.setContactForces(impulse_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  if (robot.has_floating_base()) {
    ID.head(robot.dim_passive()) -= s.u_passive;
    ID.tail(robot.dimu()) -= s.u;
  }
  else {
    ID -= s.u;
  }
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual.lq() += dtau * dID_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * dID_dv.transpose() * s.beta;
  kkt_residual.la += dtau * dID_da.transpose() * s.beta;
  if (robot.has_floating_base()) {
    kkt_residual.lu() -= dtau * s.beta.tail(robot.dimu());
    kkt_residual.lu_passive -= dtau * s.beta.head(robot.dim_passive());
    kkt_residual.lu_passive += dtau * s.nu_passive;
  }
  else {
    kkt_residual.lu() -= dtau * s.beta;
  }
  Eigen::VectorXd C = Eigen::VectorXd::Zero(impulse_status.dimf());
  if (impulse_status.hasActiveContacts()) {
    Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), impulse_status.dimf());
    robot.updateKinematics(s.q, s.v, s.a);
    robot.dRNEAPartialdFext(impulse_status, dID_df);
    kkt_residual.lf() += dtau * dID_df.transpose() * s.beta;
    Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(impulse_status.dimf(), robot.dimv());
    Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(impulse_status.dimf(), robot.dimv());
    Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(impulse_status.dimf(), robot.dimv());
    robot.updateKinematics(s.q, s.v, s.a);
    robot.computeBaumgarteResidual(impulse_status, dtau, C);
    robot.computeBaumgarteDerivatives(impulse_status, dtau, dC_dq, dC_dv, dC_da);
    kkt_residual.lq() += dtau * dC_dq.transpose() * s.mu_stack();
    kkt_residual.lv() += dtau * dC_dv.transpose() * s.mu_stack();
    kkt_residual.la += dtau * dC_da.transpose() * s.mu_stack();
  }
  double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                         + kkt_residual.lx().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + dtau * dtau * ID.squaredNorm();
  if (robot.has_floating_base()) {
    kkt_error_ref += kkt_residual.lu_passive.squaredNorm();
    kkt_error_ref += dtau * dtau * s.u_passive.squaredNorm();
  }
  if (impulse_status.hasActiveContacts()) {
    kkt_error_ref += dtau * dtau * C.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  double constraint_violation_ref = kkt_residual.Fx().lpNorm<1>() + dtau * ID.lpNorm<1>();
  if (robot.has_floating_base()) {
    constraint_violation_ref += dtau * s.u_passive.lpNorm<1>();
  }
  if (impulse_status.hasActiveContacts()) {
    constraint_violation_ref += dtau * C.lpNorm<1>();
  }
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyCost(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, empty_cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, dtau, s);
  ocp.computeKKTResidual(robot, impulse_status, t, dtau, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, impulse_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                         + kkt_residual.lx().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + cd.squaredNormContactDynamicsResidual(dtau)
                         + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  if (robot.has_floating_base()) {
    kkt_error_ref += kkt_residual.lu_passive.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testComputeKKTResidualEmptyConstraints(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost) {
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, cost, empty_constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, dtau, s);
  ocp.computeKKTResidual(robot, impulse_status, t, dtau, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, impulse_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                         + kkt_residual.lx().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + cd.squaredNormContactDynamicsResidual(dtau);
  if (robot.has_floating_base()) {
    kkt_error_ref += kkt_residual.lu_passive.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testComputeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, dtau, s);
  ocp.computeKKTResidual(robot, impulse_status, t, dtau, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, s_prev.q, s, s_next, kkt_matrix, kkt_residual);
  ContactDynamics cd(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cd.linearizeContactDynamics(robot, impulse_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.Fx().squaredNorm()
                         + kkt_residual.lx().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + cd.squaredNormContactDynamicsResidual(dtau)
                         + constraints->squaredNormPrimalAndDualResidual(constraints_data);
  if (robot.has_floating_base()) {
    kkt_error_ref += kkt_residual.lu_passive.squaredNorm();
  }
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


void SplitImpulseOCPTest::testCostAndConstraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, 
    const std::shared_ptr<ImpulseCostFunction>& cost,
    const std::shared_ptr<ImpulseConstraints>& constraints) {
  const ImpulseSplitSolution s_prev = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitSolution s_next = ImpulseSplitSolution::Random(robot, impulse_status);
  const ImpulseSplitDirection d_prev = ImpulseSplitDirection::Random(robot, impulse_status);
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  const ImpulseSplitDirection d_next = ImpulseSplitDirection::Random(robot, impulse_status);
  SplitImpulseOCP ocp(robot, cost, constraints);
  const double t = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  ocp.initConstraints(robot, 10, dtau, s);
  ocp.computeKKTResidual(robot, impulse_status, t, dtau, s_prev.q, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  const double step_size = 0.3;
  KKTMatrix kkt_matrix(robot);
  kkt_matrix.setImpulseStatus(impulse_status);
  KKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  auto cost_data = cost->createCostFunctionData(robot);
  auto constraints_data = constraints->createConstraintsData(robot, 10);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  SplitSolution s_tmp = s;
  s_tmp.integrate(robot, step_size, d);
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  double cost_ref = 0;
  cost_ref += cost->l(robot, cost_data, t, dtau, s_tmp);
  cost_ref += constraints->costSlackBarrier(constraints_data, step_size);
  constraints->computePrimalAndDualResidual(robot, constraints_data, dtau, s_tmp);
  stateequation::ComputeForwardEulerResidual(robot, step_size, dtau, s_tmp,  
                                             s_next.q, s_next.v, d_next.dq(), 
                                             d_next.dv(), kkt_residual);
  ContactDynamics cd(robot);
  cd.computeContactDynamicsResidual(robot, impulse_status, dtau, s_tmp);
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
  if (!impulse_status.hasActiveContacts()) {
    impulse_status.activateContact(0);
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