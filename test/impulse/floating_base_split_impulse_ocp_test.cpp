#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


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
    q_prev = Eigen::VectorXd::Zero(robot.dimq());
    robot.generateFeasibleConfiguration(q_prev);
    s = SplitSolution::Random(robot, contact_status);
    s_next = SplitSolution::Random(robot, contact_status);
    s_tmp = SplitSolution::Random(robot, contact_status);
    d = SplitDirection::Random(robot, contact_status);
    d_next = SplitDirection::Random(robot, contact_status);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    auto joint_cost = std::make_shared<JointSpaceCost>(robot);
    const Eigen::VectorXd q_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    Eigen::VectorXd q_ref = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_ref);
    const Eigen::VectorXd v_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd v_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd a_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd a_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd u_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd u_ref = Eigen::VectorXd::Random(robot.dimv());
    const Eigen::VectorXd qf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv()).array().abs();
    joint_cost->set_q_weight(q_weight);
    joint_cost->set_q_ref(q_ref);
    joint_cost->set_v_weight(v_weight);
    joint_cost->set_v_ref(v_ref);
    joint_cost->set_a_weight(a_weight);
    joint_cost->set_a_ref(a_ref);
    joint_cost->set_u_weight(u_weight);
    joint_cost->set_u_ref(u_ref);
    joint_cost->set_qf_weight(qf_weight);
    joint_cost->set_vf_weight(vf_weight);
    cost = std::make_shared<CostFunction>();
    cost->push_back(joint_cost);
    cost_data = cost->createCostFunctionData(robot);
    constraints = std::make_shared<Constraints>();
    auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
    auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
    auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
    auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
    constraints->push_back(joint_upper_limit); 
    constraints->push_back(joint_lower_limit);
    constraints->push_back(velocity_lower_limit); 
    constraints->push_back(velocity_upper_limit);
    constraints_data = constraints->createConstraintsData(robot, 2);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual = KKTResidual(robot);
    kkt_matrix.setContactStatus(contact_status);
    kkt_residual.setContactStatus(contact_status);
    robot_dynamics = RobotDynamics(robot);
    gain = RiccatiGain(robot);
    factorizer = RiccatiMatrixFactorizer(robot);
    inverter = RiccatiMatrixInverter(robot);
    gain.setContactStatus(contact_status);
    inverter.setContactStatus(contact_status);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  ContactStatus contact_status;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  Eigen::VectorXd q_prev;
  SplitSolution s, s_next, s_tmp;
  SplitDirection d, d_next;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  RobotDynamics robot_dynamics;
  RiccatiGain gain;
  RiccatiMatrixFactorizer factorizer;
  RiccatiMatrixInverter inverter;
};


TEST_F(FloatingBaseSplitImpulseOCPTest, isFeasible) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  EXPECT_EQ(constraints->isFeasible(robot, constraints_data, s),
            ocp.isFeasible(robot, s));
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  contact_status.setContactStatus(std::vector<bool>({false, false, false, false}));
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  s.mu_stack().setZero();
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = kkt_matrix.Fqq.transpose() * s_next.lmd 
                        + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  const double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                                + kkt_residual.Fv().squaredNorm()
                                + kkt_residual.lq().squaredNorm()
                                + kkt_residual.lv().squaredNorm()
                                + kkt_residual.la().squaredNorm()
                                + dtau*dtau*s.u.head(6).squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  const auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormStateEquationAndInverseDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  contact_status.setContactStatus(std::vector<bool>({false, false, false, false}));
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  s.mu_stack().setZero();
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = kkt_matrix.Fqq.transpose() * s_next.lmd + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  kkt_residual.lq() += dtau * du_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * du_dv.transpose() * s.beta;
  kkt_residual.la() += dtau * du_da.transpose() * s.beta;
  kkt_residual.lu -= dtau * s.beta;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm()
                         + kkt_residual.lu.squaredNorm()
                         + dtau*dtau*kkt_residual.u_res.squaredNorm()
                         + dtau*dtau*s.u.head(6).squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  SplitImpulseOCP ocp(robot, empty_cost, empty_constraints);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = kkt_matrix.Fqq.transpose() * s_next.lmd + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, du_df);
  kkt_residual.lq() += dtau * du_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * du_dv.transpose() * s.beta;
  kkt_residual.la() += dtau * du_da.transpose() * s.beta;
  kkt_residual.lf() += dtau * du_df.transpose() * s.beta;
  kkt_residual.lu -= dtau * s.beta;
  robot.computeBaumgarteResidual(contact_status, dtau, dtau, kkt_residual.C().tail(contact_status.dimf()));
  robot.computeBaumgarteDerivatives(contact_status, dtau, dtau, kkt_matrix.Cq().bottomRows(contact_status.dimf()), 
                                    kkt_matrix.Cv().bottomRows(contact_status.dimf()), 
                                    kkt_matrix.Ca().bottomRows(contact_status.dimf()));
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_stack();
  kkt_residual.lu.head(6) += dtau * s.mu_stack().head(6);
  kkt_residual.C().head(6) = dtau * s.u.head(6);
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu.squaredNorm()
                         + dtau*dtau*kkt_residual.u_res.squaredNorm()
                         + kkt_residual.C().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + dtau * kkt_residual.u_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitImpulseOCP ocp(robot, empty_cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u, kkt_residual.lu);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.linearizeRobotDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la().squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu.squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += robot_dynamics.squaredNormRobotDynamicsResidual(dtau, kkt_residual);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitImpulseOCP ocp(robot, cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.linearizeRobotDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la().squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu.squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += robot_dynamics.squaredNormRobotDynamicsResidual(dtau, kkt_residual);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, KKTErrorNorm) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  if (contact_status.hasActiveContacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, 
                                    kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u, kkt_residual.lu);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.linearizeRobotDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la().squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu.squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += robot_dynamics.squaredNormRobotDynamicsResidual(dtau, kkt_residual);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, costAndConstraintViolation) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau, dtau, kkt_residual.C_contacts());
  constraints->computePrimalAndDualResidual(robot, constraints_data, dtau, s);
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() + kkt_residual.C_contacts().lpNorm<1>()
          + dtau * s.u.head(6).lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, costAndConstraintViolationWithStepSize) {
  SplitImpulseOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double step_size = 0.3;
  const auto pair = ocp.costAndConstraintViolation(robot, contact_status, 
                                                   step_size, t, dtau, s, d, 
                                                   s_next, d_next); 
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.a = s.a + step_size * d.da();
  s_tmp.f_stack() = s.f_stack() + step_size * d.df();
  s_tmp.set_f();
  robot.setContactForces(contact_status, s_tmp.f);
  s_tmp.u = s.u + step_size * d.du;
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  const double cost_ref 
    = cost->l(robot, cost_data, t, dtau, s_tmp) 
        + constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_NEAR(pair.first, cost_ref, dtau);
  robot.subtractConfiguration(s_tmp.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s_tmp.v - step_size * d_next.dq();
  kkt_residual.Fv() = s_tmp.v + dtau * s_tmp.a - s_next.v - step_size * d_next.dv();
  robot.setContactForces(contact_status, s_tmp.f);
  robot.RNEA(s_tmp.q, s_tmp.v, s_tmp.a, kkt_residual.u_res);
  kkt_residual.u_res -= s_tmp.u;
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  robot.computeBaumgarteResidual(contact_status, dtau, dtau, kkt_residual.C_contacts());
  constraints->computePrimalAndDualResidual(robot, constraints_data, dtau, s_tmp);
  double violation_ref = 0;
  violation_ref += kkt_residual.Fx().lpNorm<1>();
  violation_ref += dtau * kkt_residual.u_res.lpNorm<1>();
  violation_ref += kkt_residual.C_contacts().lpNorm<1>();
  violation_ref += dtau * s_tmp.u.head(6).lpNorm<1>();
  violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitImpulseOCPTest, riccatiRecursion) {
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf() + robot.dim_passive();
  const int dim_passive = robot.dim_passive();
  std::cout << "dimf = " << dimf << std::endl;
  std::cout << "dimc = " << dimc << std::endl;
  s.setContactStatus(contact_status);
  SplitImpulseOCP ocp(robot, cost, constraints);
  while (!ocp.isFeasible(robot, s)) {
    s.v = Eigen::VectorXd::Random(dimv);
    s.u = Eigen::VectorXd::Random(dimv);
  }
  ASSERT_TRUE(ocp.isFeasible(robot, s));
  ocp.initConstraints(robot, 5, dtau, s);
  ocp.linearizeOCP(robot, contact_status, t, dtau, q_prev, s, s_next);
  const auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s); 
  RiccatiFactorization riccati_next(robot);
  const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
  const Eigen::MatrixXd P = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(2*dimv, 2*dimv);
  riccati_next.Pqq = P.topLeftCorner(dimv, dimv);
  riccati_next.Pqv = P.topRightCorner(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  riccati_next.Pvv = P.bottomRightCorner(dimv, dimv);
  RiccatiFactorization riccati(robot);
  ocp.backwardRiccatiRecursion(dtau, riccati_next, riccati);
  ocp.forwardRiccatiRecursion(dtau, d, d_next);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->lq(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lv(robot, cost_data, t, dtau, s, kkt_residual);
  cost->la(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lf(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  Eigen::VectorXd q_res = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd v_res = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(dimv);
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(dimv, dimf);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, kkt_matrix.Fqq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, kkt_matrix.Fqq_prev);
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.dRNEAPartialdFext(contact_status, du_df);
  robot.computeBaumgarteResidual(contact_status, dtau, dtau, 
                                 kkt_residual.C_contacts());
  robot.computeBaumgarteDerivatives(contact_status, dtau, dtau, 
                                    kkt_matrix.Cq_contacts(), 
                                    kkt_matrix.Cv_contacts(), 
                                    kkt_matrix.Ca_contacts());
  kkt_residual.C_floating_base() = dtau * (s.u.head(dim_passive)+kkt_residual.u_res.head(dim_passive));
  kkt_matrix.Cq_floating_base() = dtau * du_dq.topRows(dim_passive);
  kkt_matrix.Cv_floating_base() = dtau * du_dv.topRows(dim_passive);
  kkt_matrix.Ca_floating_base() = dtau * du_da.topRows(dim_passive);
  kkt_matrix.Cf_floating_base() = dtau * du_df.topRows(dim_passive);
  Eigen::MatrixXd Cu = Eigen::MatrixXd::Zero(dimc, dimv);
  Cu.topLeftCorner(dim_passive, dim_passive) 
      = dtau * Eigen::MatrixXd::Identity(dim_passive, dim_passive);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s.u,
                                   kkt_residual.lu);
  cost->luu(robot, cost_data, t, dtau, s.u, kkt_matrix.Quu);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s.u, 
                                    kkt_matrix.Quu, kkt_residual.lu);
  const Eigen::VectorXd lu_condensed = kkt_residual.lu + kkt_matrix.Quu * kkt_residual.u_res;
  kkt_residual.lq() += du_dq.transpose() * lu_condensed;
  kkt_residual.lv() += du_dv.transpose() * lu_condensed;
  kkt_residual.la() += du_da.transpose() * lu_condensed;
  kkt_residual.lf() += du_df.transpose() * lu_condensed;
  kkt_residual.lq() += kkt_matrix.Fqq.transpose() * s_next.lmd + kkt_matrix.Fqq_prev.transpose() * s.lmd;
  kkt_residual.lv() += dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() += dtau * s_next.gmm;
  constraints->augmentDualResidual(robot, constraints_data, dtau, s,
                                   kkt_residual);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_stack();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_stack();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_stack();
  kkt_residual.lf() += kkt_matrix.Cf().transpose() * s.mu_stack();
  kkt_residual.lu -= dtau * s.beta;
  kkt_residual.lu += Cu.transpose() * s.mu_stack();
  kkt_matrix.Qqq() = du_dq.transpose() * kkt_matrix.Quu * du_dq;
  kkt_matrix.Qqv() = du_dq.transpose() * kkt_matrix.Quu * du_dv;
  kkt_matrix.Qqa() = du_dq.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qvv() = du_dv.transpose() * kkt_matrix.Quu * du_dv;
  kkt_matrix.Qva() = du_dv.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qaa() = du_da.transpose() * kkt_matrix.Quu * du_da;
  kkt_matrix.Qqf() = du_dq.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qvf() = du_dv.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qaf() = du_da.transpose() * kkt_matrix.Quu * du_df;
  kkt_matrix.Qff() = du_df.transpose() * kkt_matrix.Quu * du_df;
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, 
                                    kkt_matrix, kkt_residual);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix);
  factorizer.setStateEquationDerivative(kkt_matrix.Fqq);
  // factorizer.setStateEquationDerivativeInverse(kkt_matrix.Fqq_prev);
  inverter.setContactStatus(contact_status);
  gain.setContactStatus(contact_status);
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();
  if (contact_status.dimf() > 0) {
    kkt_matrix.Qfq() = kkt_matrix.Qqf().transpose();
    kkt_matrix.Qfv() = kkt_matrix.Qvf().transpose();
    kkt_matrix.Qfa() = kkt_matrix.Qaf().transpose();
  }
  // std::cout << "in test" << std::endl;
  // std::cout << std::setprecision(20) << "Qqq\n" << kkt_matrix.Qqq() << std::endl;
  // std::cout << std::setprecision(20) << "Qqv\n" << kkt_matrix.Qqv() << std::endl;
  // std::cout << std::setprecision(20) << "Qqa\n" << kkt_matrix.Qqa() << std::endl;
  // std::cout << std::setprecision(20) << "Qqf\n" << kkt_matrix.Qqf() << std::endl;
  // std::cout << std::setprecision(20) << "Qvq\n" << kkt_matrix.Qvq() << std::endl;
  // std::cout << std::setprecision(20) << "Qvv\n" << kkt_matrix.Qvv() << std::endl;
  // std::cout << std::setprecision(20) << "Qva\n" << kkt_matrix.Qva() << std::endl;
  // std::cout << std::setprecision(20) << "Qvf\n" << kkt_matrix.Qvf() << std::endl;
  // std::cout << std::setprecision(20) << "Qaq\n" << kkt_matrix.Qaq() << std::endl;
  // std::cout << std::setprecision(20) << "Qav\n" << kkt_matrix.Qav() << std::endl;
  // std::cout << std::setprecision(20) << "Qaa\n" << kkt_matrix.Qaa() << std::endl;
  // std::cout << std::setprecision(20) << "Qaf\n" << kkt_matrix.Qaf() << std::endl;
  // std::cout << std::setprecision(20) << "Qfq\n" << kkt_matrix.Qfq() << std::endl;
  // std::cout << std::setprecision(20) << "Qfv\n" << kkt_matrix.Qfv() << std::endl;
  // std::cout << std::setprecision(20) << "Qfa\n" << kkt_matrix.Qfa() << std::endl;
  // std::cout << std::setprecision(20) << "Qff\n" << kkt_matrix.Qff() << std::endl;
  // std::cout << std::setprecision(20) << "Cq\n" << kkt_matrix.Cq() << std::endl;
  // std::cout << std::setprecision(20) << "Cv\n" << kkt_matrix.Cv() << std::endl;
  // std::cout << std::setprecision(20) << "Ca\n" << kkt_matrix.Ca() << std::endl;
  // std::cout << std::setprecision(20) << "Cf\n" << kkt_matrix.Cf() << std::endl;
  // std::cout << std::setprecision(20) << "lq\n" << kkt_residual.lq() << std::endl;
  // std::cout << std::setprecision(20) << "lv\n" << kkt_residual.lv() << std::endl;
  // std::cout << std::setprecision(20) << "la\n" << kkt_residual.la() << std::endl;
  // std::cout << std::setprecision(20) << "lf\n" << kkt_residual.lf() << std::endl;
  // std::cout << std::setprecision(20) << "C\n" << kkt_residual.C() << std::endl;

  factorizer.factorize_F(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                         riccati_next.Pvq, riccati_next.Pvv, 
                         kkt_matrix.Qqq(), kkt_matrix.Qqv(), 
                         kkt_matrix.Qvq(), kkt_matrix.Qvv());
  factorizer.factorize_H(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                         kkt_matrix.Qqa(), kkt_matrix.Qva());
  factorizer.factorize_G(dtau, riccati_next.Pvv, kkt_matrix.Qaa());
  factorizer.factorize_la(dtau, riccati_next.Pvq, riccati_next.Pvv, 
                          kkt_residual.Fq(), kkt_residual.Fv(), 
                          riccati_next.sv, kkt_residual.la());
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();

  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(dimv+dimf+dimc, dimv+dimf+dimc);
  G.topLeftCorner(dimv+dimf, dimv+dimf) = kkt_matrix.Qafaf();
  G.topRightCorner(dimv+dimf, dimc) = kkt_matrix.Caf().transpose();
  G.bottomLeftCorner(dimc, dimv+dimf) = kkt_matrix.Caf();
  const Eigen::MatrixXd Ginv = G.inverse();
  // Eigen::MatrixXd Ginv = Eigen::MatrixXd::Zero(dimv+dimf+dimc, dimv+dimf+dimc);
  // inverter.invert(kkt_matrix.Qafaf(), kkt_matrix.Caf(), Ginv);
  Eigen::MatrixXd Qqvaf = Eigen::MatrixXd::Zero(2*dimv, dimv+dimf);
  Qqvaf.topLeftCorner(dimv, dimv) = kkt_matrix.Qqa();
  Qqvaf.topRightCorner(dimv, dimf) = kkt_matrix.Qqf();
  Qqvaf.bottomLeftCorner(dimv, dimv) = kkt_matrix.Qva();
  Qqvaf.bottomRightCorner(dimv, dimf) = kkt_matrix.Qvf();
  EXPECT_TRUE(Qqvaf.transpose().isApprox(kkt_matrix.Qafqv()));
  gain.computeFeedbackGain(Ginv, Qqvaf.transpose(), kkt_matrix.Cqv());
  gain.computeFeedforward(Ginv, kkt_residual.laf(), kkt_residual.C());
  Eigen::MatrixXd Pqq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pqv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pvq_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::MatrixXd Pvv_ref = Eigen::MatrixXd::Zero(dimv, dimv);
  Eigen::VectorXd sq_ref = Eigen::VectorXd::Zero(dimv);
  Eigen::VectorXd sv_ref = Eigen::VectorXd::Zero(dimv);
  Pqq_ref = kkt_matrix.Qqq();
  Pqq_ref += gain.Kaq().transpose() * kkt_matrix.Qqa().transpose();
  Pqv_ref = kkt_matrix.Qqv();
  Pqv_ref += gain.Kaq().transpose() * kkt_matrix.Qva().transpose();
  Pvq_ref = kkt_matrix.Qvq();
  Pvq_ref += gain.Kav().transpose() * kkt_matrix.Qqa().transpose();
  Pvv_ref = kkt_matrix.Qvv();
  Pvv_ref += gain.Kav().transpose() * kkt_matrix.Qva().transpose();
  Pqq_ref += gain.Kfq().transpose() * kkt_matrix.Qqf().transpose();
  Pqv_ref += gain.Kfq().transpose() * kkt_matrix.Qvf().transpose();
  Pvq_ref += gain.Kfv().transpose() * kkt_matrix.Qqf().transpose();
  Pvv_ref += gain.Kfv().transpose() * kkt_matrix.Qvf().transpose();
  Pqq_ref += gain.Kmuq().transpose() * kkt_matrix.Cq();
  Pqv_ref += gain.Kmuq().transpose() * kkt_matrix.Cv();
  Pvq_ref += gain.Kmuv().transpose() * kkt_matrix.Cq();
  Pvv_ref += gain.Kmuv().transpose() * kkt_matrix.Cv();
  sq_ref = riccati_next.sq - kkt_residual.lq();
  sq_ref -= riccati_next.Pqq * kkt_residual.Fq();
  sq_ref -= riccati_next.Pqv * kkt_residual.Fv();
  sq_ref -= kkt_matrix.Qqa() * gain.ka();
  sv_ref = dtau * riccati_next.sq + riccati_next.sv - kkt_residual.lv();
  sv_ref -= dtau * riccati_next.Pqq * kkt_residual.Fq();
  sv_ref -= riccati_next.Pvq * kkt_residual.Fq();
  sv_ref -= dtau * riccati_next.Pqv * kkt_residual.Fv();
  sv_ref -= riccati_next.Pvv * kkt_residual.Fv();
  sv_ref -= kkt_matrix.Qva() * gain.ka();
  sq_ref -= kkt_matrix.Qqf() * gain.kf();
  sv_ref -= kkt_matrix.Qvf() * gain.kf();
  sq_ref -= kkt_matrix.Cq().transpose() * gain.kmu();
  sv_ref -= kkt_matrix.Cv().transpose() * gain.kmu();
  // factorizer.correctP(Pqq_ref, Pqv_ref);
  // factorizer.correct_s(sq_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(Pqq_ref, 1.0e-10));
  EXPECT_TRUE(riccati.Pqv.isApprox(Pqv_ref, 1.0e-10));
  EXPECT_TRUE(riccati.Pvq.isApprox(Pvq_ref, 1.0e-10));;
  EXPECT_TRUE(riccati.Pvv.isApprox(Pvv_ref, 1.0e-10));
  EXPECT_TRUE(riccati.sq.isApprox(sq_ref, 1.0e-10));
  EXPECT_TRUE(riccati.sv.isApprox(sv_ref, 1.0e-10));
  std::cout << riccati.Pqq - Pqq_ref << std::endl;
  std::cout << riccati.Pqv - Pqv_ref << std::endl;
  std::cout << riccati.Pvq - Pvq_ref << std::endl;
  std::cout << riccati.Pvv - Pvv_ref << std::endl;
  std::cout << riccati.sq - sq_ref << std::endl;
  std::cout << riccati.sv - sv_ref << std::endl;
  EXPECT_TRUE(riccati.Pqq.transpose().isApprox(riccati.Pqq));
  EXPECT_TRUE(riccati.Pqv.transpose().isApprox(riccati.Pvq));
  EXPECT_TRUE(riccati.Pvv.transpose().isApprox(riccati.Pvv));
  const Eigen::VectorXd da_ref = gain.ka() + gain.Kaq() * d.dq() + gain.Kav() * d.dv(); 
  const Eigen::VectorXd dq_next_ref = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  const Eigen::VectorXd dv_next_ref = d.dv() + dtau * d.da() + kkt_residual.Fv();
  EXPECT_TRUE(d_next.dq().isApprox(dq_next_ref));
  EXPECT_TRUE(d_next.dv().isApprox(dv_next_ref));
  ocp.computeCondensedDirection(robot, dtau, s, d);
  EXPECT_TRUE(d.df().isApprox(gain.kf()+gain.Kfq()*d.dq()+gain.Kfv()*d.dv(), 1.0e-10));
  EXPECT_TRUE(d.dmu().isApprox(gain.kmu()+gain.Kmuq()*d.dq()+gain.Kmuv()*d.dv(), 1.0e-10));
  Eigen::MatrixXd Kuq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Kuv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  ocp.getStateFeedbackGain(Kuq, Kuv);
  Eigen::MatrixXd Kuq_ref = du_dq + du_da * gain.Kaq();
  if (dimf > 0) {
    Kuq_ref += du_df * gain.Kfq();
  }
  Eigen::MatrixXd Kuv_ref = du_dv + du_da * gain.Kav();
  if (dimf > 0) {
    Kuv_ref += du_df * gain.Kfv();
  }
  EXPECT_TRUE(Kuq.isApprox(Kuq_ref, 1.0e-10));
  EXPECT_TRUE(Kuv.isApprox(Kuv_ref, 1.0e-10));
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const auto pair_ref = ocp.costAndConstraintViolation(robot, t, dtau, s); 
  EXPECT_DOUBLE_EQ(pair.first, pair_ref.first);
  EXPECT_DOUBLE_EQ(pair.second, pair_ref.second);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}