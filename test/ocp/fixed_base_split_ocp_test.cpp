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
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"


namespace idocp {

class FixedBaseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    std::vector<int> contact_frames = {18};
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
    contact_dynamics = ContactDynamics(robot);
    gain = RiccatiGain(robot);
    factorizer = RiccatiMatrixFactorizer(robot);
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
  ContactDynamics contact_dynamics;
  RiccatiGain gain;
  RiccatiMatrixFactorizer factorizer;
};


TEST_F(FixedBaseSplitOCPTest, isFeasible) {
  SplitOCP ocp(robot, cost, constraints);
  EXPECT_EQ(constraints->isFeasible(robot, constraints_data, s),
            ocp.isFeasible(robot, s));
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  contact_status.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la = dtau * s_next.gmm;
  const double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                                + kkt_residual.Fv().squaredNorm()
                                + kkt_residual.lq().squaredNorm()
                                + kkt_residual.lv().squaredNorm()
                                + kkt_residual.la.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  const auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormStateEquationAndInverseDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  contact_status.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la = dtau * s_next.gmm;
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  kkt_residual.lq() += dtau * dID_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * dID_dv.transpose() * s.beta;
  kkt_residual.la += dtau * dID_da.transpose() * s.beta;
  kkt_residual.lu() -= dtau * s.beta;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + dtau * dtau * ID.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormStateEquationAndContactDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  contact_status.setContactStatus(std::vector<bool>({true}));
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.Fq() = s.q + dtau * s.v - s_next.q;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la = dtau * s_next.gmm;
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dID_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dID_df = Eigen::MatrixXd::Zero(robot.dimv(), contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, dID_dq, dID_dv, dID_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(contact_status, dID_df);
  kkt_residual.lq() += dtau * dID_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * dID_dv.transpose() * s.beta;
  kkt_residual.la += dtau * dID_da.transpose() * s.beta;
  kkt_residual.lf() += dtau * dID_df.transpose() * s.beta;
  kkt_residual.lu() -= dtau * s.beta;
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  Eigen::MatrixXd dC_dq = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_dv = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  Eigen::MatrixXd dC_da = Eigen::MatrixXd::Zero(contact_status.dimf(), robot.dimv());
  robot.computeBaumgarteResidual(contact_status, dtau, C);
  robot.computeBaumgarteDerivatives(contact_status, dtau, dC_dq, dC_dv, dC_da);
  kkt_residual.lq() += dtau * dC_dq.transpose() * s.mu_stack();
  kkt_residual.lv() += dtau * dC_dv.transpose() * s.mu_stack();
  kkt_residual.la += dtau * dC_da.transpose() * s.mu_stack();
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la.squaredNorm()
                         + kkt_residual.lf().squaredNorm()
                         + kkt_residual.lu().squaredNorm()
                         + dtau * dtau * ID.squaredNorm()
                         + dtau * dtau * C.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + dtau * ID.lpNorm<1>()
                                + dtau * C.lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitOCP ocp(robot, empty_cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  contact_dynamics.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += contact_dynamics.squaredNormContactDynamicsResidual(dtau);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitOCP ocp(robot, cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  ocp.computeKKTResidual(robot, contact_status, t, dtau, q_prev, s, s_next);
  const double kkt_error = ocp.squaredNormKKTResidual(dtau);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setContactStatus(contact_status);
  s.setContactStatus(contact_status);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  contact_dynamics.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += contact_dynamics.squaredNormContactDynamicsResidual(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, KKTErrorNorm) {
  SplitOCP ocp(robot, cost, constraints);
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
  constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
  contact_dynamics.linearizeContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = 0;
  kkt_error_ref += kkt_residual.lq().squaredNorm();
  kkt_error_ref += kkt_residual.lv().squaredNorm();
  kkt_error_ref += kkt_residual.la.squaredNorm();
  kkt_error_ref += kkt_residual.lf().squaredNorm();
  kkt_error_ref += kkt_residual.lu().squaredNorm();
  kkt_error_ref += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  kkt_error_ref += contact_dynamics.squaredNormContactDynamicsResidual(dtau);
  kkt_error_ref += constraints->squaredNormPrimalAndDualResidual(constraints_data);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseSplitOCPTest, costAndConstraintViolation) {
  SplitOCP ocp(robot, cost, constraints);
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
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.setContactForces(contact_status, s.f);
  robot.RNEA(s.q, s.v, s.a, ID);
  ID -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(contact_status, dtau, C);
  constraints->computePrimalAndDualResidual(robot, constraints_data, dtau, s);
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * ID.lpNorm<1>() + dtau * C.lpNorm<1>()
          + constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseSplitOCPTest, costAndConstraintViolationWithStepSize) {
  SplitOCP ocp(robot, cost, constraints);
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
  s_tmp.u = s.u + step_size * d.du();
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s_tmp) 
          + constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s_tmp.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s_tmp.v - step_size * d_next.dq();
  kkt_residual.Fv() = s_tmp.v + dtau * s_tmp.a - s_next.v - step_size * d_next.dv();
  Eigen::VectorXd ID = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::VectorXd C = Eigen::VectorXd::Zero(contact_status.dimf());
  robot.setContactForces(contact_status, s_tmp.f);
  robot.RNEA(s_tmp.q, s_tmp.v, s_tmp.a, ID);
  ID -= s_tmp.u;
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  robot.computeBaumgarteResidual(contact_status, dtau, C);
  constraints->computePrimalAndDualResidual(robot, constraints_data, dtau, s_tmp);
  double violation_ref = 0;
  violation_ref += kkt_residual.Fx().lpNorm<1>();
  violation_ref += dtau * ID.lpNorm<1>();
  violation_ref += dtau * C.lpNorm<1>();
  violation_ref += constraints->l1NormPrimalResidual(constraints_data);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


// TEST_F(FixedBaseSplitOCPTest, riccatiRecursion) {
//   const int dimv = robot.dimv();
//   RiccatiFactorization riccati_next(robot);
//   const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(2*dimv, 2*dimv);
//   const Eigen::MatrixXd P = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(2*dimv, 2*dimv);
//   riccati_next.Pqq = P.topLeftCorner(dimv, dimv);
//   riccati_next.Pqv = P.topRightCorner(dimv, dimv);
//   riccati_next.Pvq = riccati_next.Pqv.transpose();
//   riccati_next.Pvv = P.bottomRightCorner(dimv, dimv);
//   s.setContactStatus(contact_status);
//   SplitOCP ocp(robot, cost, constraints);
//   ASSERT_FALSE(robot.has_floating_base());
//   ASSERT_TRUE(robot.dim_passive() == 0);
//   while (!ocp.isFeasible(robot, s)) {
//     s.v = Eigen::VectorXd::Random(dimv);
//     s.u = Eigen::VectorXd::Random(dimv);
//   }
//   ASSERT_TRUE(ocp.isFeasible(robot, s));
//   ocp.initConstraints(robot, 5, dtau, s);
//   ocp.linearizeOCP(robot, contact_status, t, dtau, q_prev, s, s_next);
//   RiccatiFactorization riccati(robot);
//   ocp.backwardRiccatiRecursion(dtau, riccati_next, riccati);
//   ocp.forwardRiccatiRecursion(dtau, d, d_next);
//   const auto pair = ocp.costAndConstraintViolation(robot, t, dtau, s); 
//   robot.updateKinematics(s.q, s.v, s.a);
//   cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
//   cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix);
//   constraints->augmentDualResidual(robot, constraints_data, dtau, s, kkt_residual);
//   constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, kkt_matrix, kkt_residual);
//   stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, kkt_matrix, kkt_residual);
//   contact_dynamics.condenseContactDynamics(robot, contact_status, dtau, s, kkt_matrix, kkt_residual);
//   factorizer.factorizeMatrices(riccati_next, dtau, kkt_matrix, kkt_residual);
//   gain.computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual);
//   RiccatiFactorization riccati_ref(robot);
//   factorizer.factorizeRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, gain, riccati_ref);
//   EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
//   EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
//   EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
//   EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
//   EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
//   EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
//   EXPECT_TRUE(riccati.Pqq.transpose().isApprox(riccati.Pqq));
//   EXPECT_TRUE(riccati.Pqv.transpose().isApprox(riccati.Pvq));
//   EXPECT_TRUE(riccati.Pvv.transpose().isApprox(riccati.Pvv));
//   SplitDirection d_ref = SplitDirection(robot);
//   d_ref.setContactStatus(contact_status);
//   const Eigen::VectorXd du_ref = gain.K * d.dx() + gain.k;
//   const Eigen::VectorXd dx_next_ref = kkt_matrix.Fxx() * d.dx() + kkt_matrix.Fxu() * d.du() + kkt_residual.Fx();
//   EXPECT_TRUE(du_ref.isApprox(d.du()));
//   EXPECT_TRUE(dx_next_ref.isApprox(d_next.dx()));
// }

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}