#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/impulse_constraints.hpp"

#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/impulse/impulse_riccati_matrix_factorizer.hpp"


namespace idocp {

class FloatingBaseSplitImpulseOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/anymal/anymal.urdf";
    std::vector<int> contact_frames = {14, 24, 34, 44};
    robot = Robot(urdf, contact_frames);
    std::random_device rnd;
    std::vector<bool> is_contact_active;
    for (const auto frame : contact_frames) {
      is_contact_active.push_back(rnd()%2==0);
    }
    contact_status = ContactStatus(robot.max_point_contacts());
    contact_status.setContactStatus(is_contact_active);
    if (!contact_status.hasActiveContacts()) {
      contact_status.activateContact(0);
    }
    q_prev = Eigen::VectorXd::Zero(robot.dimq());
    robot.generateFeasibleConfiguration(q_prev);
    s = ImpulseSplitSolution::Random(robot, contact_status);
    s_next = SplitSolution::Random(robot, contact_status);
    s_tmp = ImpulseSplitSolution::Random(robot, contact_status);
    d = ImpulseSplitDirection::Random(robot, contact_status);
    d_next = SplitDirection::Random(robot, contact_status);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    auto joint_cost = std::make_shared<JointSpaceImpulseCost>(robot);
    auto impulse_cost = std::make_shared<ImpulseForceCost>(robot);
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
    std::vector<Eigen::Vector3d> f_weight, f_ref;
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      f_weight.push_back(Eigen::Vector3d::Constant(0.01));
      f_ref.push_back(Eigen::Vector3d::Zero());
    }
    impulse_cost->set_f_weight(f_weight);
    impulse_cost->set_f_ref(f_ref);
    cost = std::make_shared<ImpulseCostFunction>();
    cost->push_back(joint_cost);
    cost->push_back(impulse_cost);
    cost_data = cost->createCostFunctionData(robot);
    constraints = std::make_shared<ImpulseConstraints>();
    constraints_data = constraints->createConstraintsData(robot);
    kkt_matrix = ImpulseKKTMatrix(robot);
    kkt_residual = ImpulseKKTResidual(robot);
    kkt_matrix.setContactStatus(contact_status);
    kkt_residual.setContactStatus(contact_status);
    impulse_dynamics = ImpulseDynamicsForwardEuler(robot);
    factorizer = ImpulseRiccatiMatrixFactorizer(robot);
  }

  virtual void TearDown() {
  }

  double t;
  std::string urdf;
  Robot robot;
  ContactStatus contact_status;
  std::shared_ptr<ImpulseCostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<ImpulseConstraints> constraints;
  ConstraintsData constraints_data;
  Eigen::VectorXd q_prev;
  ImpulseSplitSolution s, s_tmp;
  SplitSolution s_next;
  ImpulseSplitDirection d;
  SplitDirection d_next;
  ImpulseKKTMatrix kkt_matrix;
  ImpulseKKTResidual kkt_residual;
  ImpulseDynamicsForwardEuler impulse_dynamics;
  ImpulseRiccatiMatrixFactorizer factorizer;
};


TEST_F(FloatingBaseSplitImpulseOCPTest, isFeasible) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  EXPECT_EQ(constraints->isFeasible(robot, constraints_data, s),
            ocp.isFeasible(robot, s));
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormStateEquation) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, s);
  s.beta.setZero();
  s.mu_stack().setZero();
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  kkt_residual.Fv() = s.v + s.dv - s_next.v;
  kkt_residual.lq() = kkt_matrix.Fqq.transpose() * s_next.lmd 
                        + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() = s_next.gmm - s.gmm;
  kkt_residual.ldv = s_next.gmm;
  Eigen::VectorXd dv_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd ddv_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.updateKinematics(s.q, s.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.ldv.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.dv_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + kkt_residual.dv_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  kkt_residual.Fv() = s.v + s.dv - s_next.v;
  kkt_residual.lq() = kkt_matrix.Fqq.transpose() * s_next.lmd 
                        + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() = s_next.gmm - s.gmm;
  kkt_residual.ldv = s_next.gmm;
  Eigen::VectorXd dv_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd ddv_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_ddv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd ddv_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.RNEAImpulseDerivatives(s.q, s.dv, ddv_dq, ddv_ddv);
  robot.updateKinematics(s.q, s.v);
  robot.dRNEAPartialdFext(contact_status, ddv_df);
  kkt_residual.lq() += ddv_dq.transpose() * s.beta;
  kkt_residual.ldv += ddv_ddv.transpose() * s.beta;
  kkt_residual.lf() += ddv_df.transpose() * s.beta;
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  robot.computeContactVelocityDerivatives(contact_status, kkt_matrix.Cq(), 
                                          kkt_matrix.Cv());
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.ldv += kkt_matrix.Cv().transpose() * s.mu_stack();
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.ldv.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.dv_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + kkt_residual.dv_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<ImpulseCostFunction>();
  SplitImpulseOCP ocp(robot, empty_cost, constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<ImpulseConstraints>();
  SplitImpulseOCP ocp(robot, cost, empty_constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNorm) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual();
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, kkt_matrix, kkt_residual);
  impulse_dynamics.linearizeImpulseDynamics(robot, contact_status, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.ldv.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += impulse_dynamics.squaredNormImpulseDynamicsResidual(kkt_residual);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, costAndConstraintViolation) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  ocp.computeKKTResidual(robot, contact_status, t, q_prev, s, s_next);
  const auto pair = ocp.costAndConstraintViolation(robot, t, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, s) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fv() = s.v + s.dv - s_next.v;
  robot.setContactForces(contact_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, kkt_residual.dv_res);
  robot.updateKinematics(s.q, s.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  constraints->computePrimalAndDualResidual(robot, constraints_data, s);
  const double violation_ref 
      = kkt_residual.Fx().lpNorm<1>()  
          + kkt_residual.dv_res.lpNorm<1>() 
          + kkt_residual.C().lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, costAndConstraintViolationWithStepSize) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, s);
  constraints->setSlackAndDual(robot, constraints_data, s);
  const double step_size = 0.3;
  const auto pair = ocp.costAndConstraintViolation(robot, contact_status, 
                                                   step_size, t, s, d, 
                                                   s_next, d_next); 
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.dv = s.dv + step_size * d.ddv;
  s_tmp.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp.set_f();
  robot.updateKinematics(s_tmp.q, s_tmp.v);
  const double cost_ref 
      = cost->l(robot, cost_data, t, s_tmp) 
          + constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s_tmp.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() -= step_size * d_next.dq();
  kkt_residual.Fv() = s_tmp.v + s_tmp.dv - s_next.v - step_size * d_next.dv();
  robot.setContactForces(contact_status, s_tmp.f);
  robot.RNEAImpulse(s_tmp.q, s_tmp.dv, kkt_residual.dv_res);
  robot.updateKinematics(s_tmp.q, s_tmp.v);
  robot.computeContactVelocityResidual(contact_status, kkt_residual.C());
  constraints->computePrimalAndDualResidual(robot, constraints_data, s_tmp);
  double violation_ref = 0;
  violation_ref += kkt_residual.Fx().lpNorm<1>();
  violation_ref += kkt_residual.dv_res.lpNorm<1>();
  violation_ref += kkt_residual.C().lpNorm<1>();
  violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, riccatiRecursion) {
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf();
  SplitImpulseOCP ocp(robot, cost, constraints);
  while (!ocp.isFeasible(robot, s)) {
    s.v = Eigen::VectorXd::Random(dimv);
    s.dv = Eigen::VectorXd::Random(dimv);
  }
  ASSERT_TRUE(ocp.isFeasible(robot, s));
  ocp.initConstraints(robot, s);
  ocp.linearizeOCP(robot, contact_status, t, q_prev, s, s_next);
  const auto pair = ocp.costAndConstraintViolation(robot, t, s); 
  RiccatiFactorization riccati_next(robot);
  const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  const Eigen::MatrixXd P = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(2*dimv, 2*dimv);
  riccati_next.Pqq = P.topLeftCorner(dimv, dimv);
  riccati_next.Pqv = P.topRightCorner(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  riccati_next.Pvv = P.bottomRightCorner(dimv, dimv);
  RiccatiFactorization riccati(robot);
  ocp.backwardRiccatiRecursion(riccati_next, riccati);
  ocp.forwardRiccatiRecursion(d, d_next);
  robot.updateKinematics(s.q, s.v);
  cost->computeStageCostHessian(robot, cost_data, t, s, kkt_matrix);
  cost->computeStageCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->setSlackAndDual(robot, constraints_data, s);
  constraints->augmentDualResidual(robot, constraints_data, s, kkt_residual);
  constraints->condenseSlackAndDual(robot, constraints_data, s, kkt_matrix, kkt_residual);
  Eigen::MatrixXd ddv_dq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd ddv_ddv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd ddv_df = Eigen::MatrixXd::Zero(dimv, dimf);
  stateequation::LinearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  impulse_dynamics.condenseImpulseDynamics(robot, contact_status, s, kkt_matrix, 
                                           kkt_residual);
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  RiccatiFactorization riccati_ref(robot);
  factorizer.factorize(kkt_matrix, kkt_residual, riccati_next, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.transpose().isApprox(riccati.Pqq));
  EXPECT_TRUE(riccati.Pqv.transpose().isApprox(riccati.Pvq));
  EXPECT_TRUE(riccati.Pvv.transpose().isApprox(riccati.Pvv));
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix.Fqq;
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq;
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv;
  const Eigen::VectorXd dx_next_ref = A * d.dx() + kkt_residual.Fx();
  EXPECT_TRUE(d_next.dq().isApprox(dx_next_ref.head(dimv)));
  EXPECT_TRUE(d_next.dv().isApprox(dx_next_ref.tail(dimv)));
  ImpulseSplitDirection d_ref = d;
  ocp.computeCondensedDirection(robot, s, d_next, d);
  impulse_dynamics.computeCondensedDirection(kkt_matrix, kkt_residual, d_next, d_ref);
  EXPECT_TRUE(d.ddv.isApprox(d_ref.ddv));
  EXPECT_TRUE(d.df().isApprox(d_ref.df()));
  EXPECT_TRUE(d.dbeta.isApprox(d_ref.dbeta));
  EXPECT_TRUE(d.dmu().isApprox(d_ref.dmu()));
  const double step_size = 0.3;
  const Eigen::VectorXd lmd_ref 
      = s.lmd + step_size * (riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq);
  const Eigen::VectorXd gmm_ref 
      = s.gmm + step_size * (riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv);
  Eigen::VectorXd q_ref = s.q;
  robot.integrateConfiguration(d.dq(), step_size, q_ref);
  const Eigen::VectorXd v_ref = s.v + step_size * d.dv();
  const Eigen::VectorXd dv_ref = s.dv + step_size * d.ddv;
  const Eigen::VectorXd f_ref = s.f_stack() + step_size * d.df();
  const Eigen::VectorXd mu_ref = s.mu_stack() + step_size * d.dmu();
  const Eigen::VectorXd beta_ref = s.beta + step_size * d.dbeta;
  ocp.updatePrimal(robot, step_size, riccati, d, s);
  EXPECT_TRUE(lmd_ref.isApprox(s.lmd));
  EXPECT_TRUE(gmm_ref.isApprox(s.gmm));
  EXPECT_TRUE(q_ref.isApprox(s.q));
  EXPECT_TRUE(v_ref.isApprox(s.v));
  EXPECT_TRUE(dv_ref.isApprox(s.dv));
  EXPECT_TRUE(f_ref.isApprox(s.f_stack()));
  EXPECT_TRUE(mu_ref.isApprox(s.mu_stack()));
  EXPECT_TRUE(beta_ref.isApprox(s.beta));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}