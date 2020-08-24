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
  RiccatiGain riccati_gain;
  RiccatiMatrixFactorizer riccati_factorizer;
  RiccatiMatrixInverter riccati_inverter;
  Eigen::MatrixXd Ginv;
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


TEST_F(FloatingBaseSplitOCPTest, linearizeOCP) {
  s.setContactStatus(robot);
  SplitOCP ocp(robot, cost, constraints);
  ocp.linearizeOCP(robot, t, dtau, s, s_next);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}