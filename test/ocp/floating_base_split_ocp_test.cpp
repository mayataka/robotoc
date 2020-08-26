#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/cost_function_data.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"

#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_matrix_factorizer.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class FloatingBaseSplitOCPTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/anymal/anymal.urdf";
    std::vector<int> contact_frames = {14, 24, 34, 44};
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, baum_a, baum_b);
    std::random_device rnd;
    for (const auto frame : contact_frames) {
      contact_status.push_back(rnd()%2==0);
    }
    robot.setContactStatus(contact_status);
    s = SplitSolution(robot);
    s.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s.q);
    s.v = Eigen::VectorXd::Random(robot.dimv());
    s.a = Eigen::VectorXd::Random(robot.dimv());
    s.f = Eigen::VectorXd::Random(robot.max_dimf());
    s.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s.lmd = Eigen::VectorXd::Random(robot.dimv());
    s.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_next = SplitSolution(robot);
    s_next.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_next.q);
    s_next.v = Eigen::VectorXd::Random(robot.dimv());
    s_next.a = Eigen::VectorXd::Random(robot.dimv());
    s_next.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_next.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_next.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_next.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_tmp= SplitSolution(robot);
    s_tmp.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_tmp.q);
    s_tmp.v = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.a = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_tmp.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_tmp.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_tmp.gmm = Eigen::VectorXd::Random(robot.dimv());
    d = SplitDirection(robot);
    d.dq() = Eigen::VectorXd::Random(robot.dimv());
    d.dv() = Eigen::VectorXd::Random(robot.dimv());
    d.da() = Eigen::VectorXd::Random(robot.dimv());
    d.df() = Eigen::VectorXd::Random(robot.dimf());
    d.du = Eigen::VectorXd::Random(robot.dimv());
    d_next = SplitDirection(robot);
    d_next.dq() = Eigen::VectorXd::Random(robot.dimv());
    d_next.dv() = Eigen::VectorXd::Random(robot.dimv());
    d_next.da() = Eigen::VectorXd::Random(robot.dimv());
    d_next.df() = Eigen::VectorXd::Random(robot.dimf());
    d_next.du = Eigen::VectorXd::Random(robot.dimv());
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    auto joint_cost = std::make_shared<JointSpaceCost>(robot);
    auto contact_cost = std::make_shared<ContactCost>(robot);
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
    const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot.max_dimf()).array().abs();
    const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot.max_dimf());
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
    contact_cost->set_f_weight(f_weight);
    contact_cost->set_f_ref(f_ref);
    cost = std::make_shared<CostFunction>();
    cost->push_back(joint_cost);
    cost->push_back(contact_cost);
    cost_data = CostFunctionData(robot);
    constraints = std::make_shared<Constraints>();
    auto joint_lower_limit = std::make_shared<JointPositionLowerLimit>(robot);
    auto joint_upper_limit = std::make_shared<JointPositionUpperLimit>(robot);
    auto velocity_lower_limit = std::make_shared<JointVelocityLowerLimit>(robot);
    auto velocity_upper_limit = std::make_shared<JointVelocityUpperLimit>(robot);
    constraints->push_back(joint_upper_limit); 
    constraints->push_back(joint_lower_limit);
    constraints->push_back(velocity_lower_limit); 
    constraints->push_back(velocity_upper_limit);
    constraints_data = constraints->createConstraintsData(robot);
    kkt_matrix = KKTMatrix(robot);
    kkt_residual = KKTResidual(robot);
    state_equation = StateEquation(robot);
    robot_dynamics = RobotDynamics(robot);
    gain = RiccatiGain(robot);
    factorizer = RiccatiMatrixFactorizer(robot);
    inverter = RiccatiMatrixInverter(robot);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::vector<bool> contact_status;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  SplitSolution s, s_next, s_tmp;
  SplitDirection d, d_next;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  StateEquation state_equation;
  RobotDynamics robot_dynamics;
  RiccatiGain gain;
  RiccatiMatrixFactorizer factorizer;
  RiccatiMatrixInverter inverter;
};


TEST_F(FloatingBaseSplitOCPTest, isFeasible) {
  s.setContactStatus(robot);
  SplitOCP ocp(robot, cost, constraints);
  EXPECT_EQ(constraints->isFeasible(robot, constraints_data, s),
            ocp.isFeasible(robot, s));
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false, false, false, false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  s.mu.setZero();
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm()
                         + dtau*dtau*s.u.head(6).squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNormStateEquationAndInverseDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false, false, false, false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  s.mu.setZero();
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
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
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitOCP ocp(robot, empty_cost, empty_constraints);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  kkt_residual.lq() = s_next.lmd - s.lmd;
  kkt_residual.lv() = dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() = dtau * s_next.gmm;
  Eigen::VectorXd u_res = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd du_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_dv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_da = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd du_df = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimf());
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.updateKinematics(s.q, s.v, s.a);
  robot.dRNEAPartialdFext(du_df);
  kkt_residual.lq() += dtau * du_dq.transpose() * s.beta;
  kkt_residual.lv() += dtau * du_dv.transpose() * s.beta;
  kkt_residual.la() += dtau * du_da.transpose() * s.beta;
  kkt_residual.lf() += dtau * du_df.transpose() * s.beta;
  kkt_residual.lu -= dtau * s.beta;
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), 
                                    kkt_matrix.Ca());
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu_active();
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu_active();
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu_active();
  kkt_residual.lu.head(6) += dtau * s.mu_active().tail(6);
  kkt_residual.C().tail(6) = dtau * s.u.head(6);
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
  auto pair = ocp.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + dtau * kkt_residual.u_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitOCP ocp(robot, empty_cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitOCP ocp(robot, cost, empty_constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitOCPTest, KKTErrorNorm) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  kkt_matrix.setContactStatus(robot);
  kkt_residual.setContactStatus(robot);
  kkt_residual.setZeroMinimum();
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, 
                                     kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual.lu);
  state_equation.linearizeForwardEuler(robot, dtau, s, s_next, kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitOCPTest, costAndViolation) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = ocp.squaredKKTErrorNorm(robot, t, dtau, s, s_next);
  const auto pair = ocp.costAndViolation(robot, t, dtau, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s) 
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  robot.subtractConfiguration(s.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s.v;
  kkt_residual.Fv() = s.v + dtau * s.a - s_next.v;
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + dtau * s.u.head(6).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitOCPTest, costAndViolationWithStepSize) {
  SplitOCP ocp(robot, cost, constraints);
  ocp.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double step_size = 0.3;
  const auto pair = ocp.costAndViolation(robot, step_size, t, dtau, s, d, s_next, d_next); 
  s_tmp.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp.q);
  s_tmp.v = s.v + step_size * d.dv();
  s_tmp.u = s.u + step_size * d.du;
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  }
  double cost_ref = 0;
  cost_ref += cost->l(robot, cost_data, t, dtau, s_tmp);
  cost_ref += constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_NEAR(pair.first, cost_ref, dtau);
  robot.subtractConfiguration(s_tmp.q, s_next.q, kkt_residual.Fq());
  kkt_residual.Fq() += dtau * s_tmp.v - step_size * d_next.dq();
  kkt_residual.Fv() = s_tmp.v + dtau * s_tmp.a - s_next.v - step_size * d_next.dv();
  robot.RNEA(s_tmp.q, s_tmp.v, s_tmp.a, kkt_residual.u_res);
  kkt_residual.u_res -= s_tmp.u;
  robot.updateKinematics(s_tmp.q, s_tmp.v, s_tmp.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + dtau * s_tmp.u.head(6).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s_tmp);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FloatingBaseSplitOCPTest, riccatiRecursion) {
  const int dimv = robot.dimv();
  const int dimf = robot.dimf();
  const int dimc = robot.dimf() + robot.dim_passive();
  const int dim_passive = robot.dim_passive();
  s.setContactStatus(robot);
  SplitOCP ocp(robot, cost, constraints);
  while (!ocp.isFeasible(robot, s)) {
    s.v = Eigen::VectorXd::Random(dimv);
    s.u = Eigen::VectorXd::Random(dimv);
  }
  ASSERT_TRUE(ocp.isFeasible(robot, s));
  ocp.initConstraints(robot, 5, dtau, s);
  ocp.linearizeOCP(robot, t, dtau, s, s_next);
  RiccatiFactorization riccati_next(robot);
  Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pqq = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pqv = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pvq = seed_mat * seed_mat.transpose();
  seed_mat = Eigen::MatrixXd::Random(robot.dimv(), robot.dimv());
  riccati_next.Pvv = seed_mat * seed_mat.transpose();
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
  robot.dIntegrateConfiguration(s.q, s.v, dtau, kkt_matrix.Fqq, kkt_matrix.Fqv);
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res.noalias() -= s.u;
  robot.RNEADerivatives(s.q, s.v, s.a, du_dq, du_dv, du_da);
  robot.dRNEAPartialdFext(du_df);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  robot.computeBaumgarteDerivatives(dtau, kkt_matrix.Cq(), kkt_matrix.Cv(), kkt_matrix.Ca());
  kkt_residual.C().tail(dim_passive) = dtau * (s.u.head(dim_passive)+kkt_residual.u_res.head(dim_passive));
  kkt_matrix.Cq().bottomRows(dim_passive) = dtau * du_dq.topRows(dim_passive);
  kkt_matrix.Cv().bottomRows(dim_passive) = dtau * du_dv.topRows(dim_passive);
  kkt_matrix.Ca().bottomRows(dim_passive) = dtau * du_da.topRows(dim_passive);
  kkt_matrix.Cf().bottomRows(dim_passive) = dtau * du_df.topRows(dim_passive);
  Eigen::MatrixXd Cu = Eigen::MatrixXd::Zero(dimc, dimv);
  Cu.bottomLeftCorner(dim_passive, dim_passive) 
      = dtau * Eigen::MatrixXd::Identity(dim_passive, dim_passive);
  kkt_residual.lu += Cu.transpose() * s.mu_active();
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  constraints->augmentDualResidual(robot, constraints_data, dtau, 
                                    kkt_residual.lu);
  cost->luu(robot, cost_data, t, dtau, s.u, kkt_matrix.Quu);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s.u, 
                                    kkt_matrix.Quu, kkt_residual.lu);
  const Eigen::VectorXd lu_condensed = kkt_residual.lu + kkt_matrix.Quu * kkt_residual.u_res;
  kkt_residual.lq() += du_dq.transpose() * lu_condensed;
  kkt_residual.lv() += du_dv.transpose() * lu_condensed;;
  kkt_residual.la() += du_da.transpose() * lu_condensed;;
  kkt_residual.lf() += du_df.transpose() * lu_condensed;
  kkt_residual.lq() += s_next.lmd - s.lmd;
  kkt_residual.lv() += dtau * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual.la() += dtau * s_next.gmm;
  constraints->augmentDualResidual(robot, constraints_data, dtau, 
                                   kkt_residual);
  kkt_residual.lq() += kkt_matrix.Cq().transpose() * s.mu.head(dimc);
  kkt_residual.lv() += kkt_matrix.Cv().transpose() * s.mu.head(dimc);
  kkt_residual.la() += kkt_matrix.Ca().transpose() * s.mu.head(dimc);
  kkt_residual.lf() += kkt_matrix.Cf().transpose() * s.mu.head(dimc);
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
  factorizer.setIntegrationSensitivities(kkt_matrix.Fqq, kkt_matrix.Fqv);
  inverter.setContactStatus(robot);
  gain.setContactStatus(robot);
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();
  if (robot.dimf() > 0) {
    kkt_matrix.Qfq() = kkt_matrix.Qqf().transpose();
    kkt_matrix.Qfv() = kkt_matrix.Qvf().transpose();
    kkt_matrix.Qfa() = kkt_matrix.Qaf().transpose();
  }
  factorizer.factorizeF(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                        riccati_next.Pvq, riccati_next.Pvv, 
                        kkt_matrix.Qqq(), kkt_matrix.Qqv(), 
                        kkt_matrix.Qvq(), kkt_matrix.Qvv());
  factorizer.factorizeH(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                        kkt_matrix.Qqa(), kkt_matrix.Qva());
  factorizer.factorizeG(dtau, riccati_next.Pvv, kkt_matrix.Qaa());
  kkt_matrix.Qaq() = kkt_matrix.Qqa().transpose();
  kkt_matrix.Qav() = kkt_matrix.Qva().transpose();
  kkt_residual.la() += dtau * riccati_next.Pvq * kkt_residual.Fq();
  kkt_residual.la() += dtau * riccati_next.Pvv * kkt_residual.Fv();
  kkt_residual.la() -= dtau * riccati_next.sv;
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
  EXPECT_TRUE(riccati.Pqq.isApprox(Pqq_ref));
  EXPECT_TRUE(riccati.Pqv.isApprox(Pqv_ref));
  EXPECT_TRUE(riccati.Pvq.isApprox(Pvq_ref));;
  EXPECT_TRUE(riccati.Pvv.isApprox(Pvv_ref));
  EXPECT_TRUE(riccati.sq.isApprox(sq_ref));
  EXPECT_TRUE(riccati.sv.isApprox(sv_ref));
  const Eigen::VectorXd da_ref = gain.ka() + gain.Kaq() * d.dq() + gain.Kav() * d.dv(); 
  const Eigen::VectorXd dq_next_ref = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  const Eigen::VectorXd dv_next_ref = d.dv() + dtau * d.da() + kkt_residual.Fv();
  EXPECT_TRUE(d_next.dq().isApprox(dq_next_ref));
  EXPECT_TRUE(d_next.dv().isApprox(dv_next_ref));
  std::cout << "C" << std::endl;
  std::cout << kkt_residual.C().transpose() << std::endl;
  std::cout << "dtau * u.head<6>" << std::endl;
  std::cout << dtau * s.u.head(6).transpose() << std::endl;
  std::cout << "baumgarte" << std::endl;
  kkt_residual.C().setZero();
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  std::cout << kkt_residual.C().head(robot.dimf()).transpose() << std::endl;
  std::cout << dtau * s.u.head(6).transpose() << std::endl;
  std::cout << "Cq" << std::endl;
  std::cout << kkt_matrix.Cq().transpose() << std::endl;
  std::cout << "Cv" << std::endl;
  std::cout << kkt_matrix.Cv().transpose() << std::endl;
  std::cout << "Cqv" << std::endl;
  std::cout << kkt_matrix.Cqv().transpose() << std::endl;
  std::cout << "Ca" << std::endl;
  std::cout << kkt_matrix.Ca().transpose() << std::endl;
  std::cout << "Cf" << std::endl;
  std::cout << kkt_matrix.Cf().transpose() << std::endl;
  std::cout << "Caf" << std::endl;
  std::cout << kkt_matrix.Caf().transpose() << std::endl;

  ocp.computeCondensedDirection(robot, dtau, d);
  EXPECT_TRUE(d.df().isApprox(gain.kf()+gain.Kfq()*d.dq()+gain.Kfv()*d.dv()));
  EXPECT_TRUE(d.dmu().isApprox(gain.kmu()+gain.Kmuq()*d.dq()+gain.Kmuv()*d.dv()));

  std::cout << "Qafaf" << std::endl;
  std::cout << kkt_matrix.Qafaf() << std::endl;
  std::cout << "Qafqv" << std::endl;
  std::cout << kkt_matrix.Qafqv() << std::endl;
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}