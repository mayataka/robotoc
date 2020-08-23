#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
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

class FixedBaseTerminalParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    std::vector<int> contact_frames = {18};
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, baum_a, baum_b);
    std::random_device rnd;
    contact_status.push_back(rnd()%2==0);
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
    s_tmp = SplitSolution(robot);
    s_old = SplitSolution(robot);
    s_old.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_old.q);
    s_old.v = Eigen::VectorXd::Random(robot.dimv());
    s_old.a = Eigen::VectorXd::Random(robot.dimv());
    s_old.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_old.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_old.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_old.gmm = Eigen::VectorXd::Random(robot.dimv());
    s_new = SplitSolution(robot);
    s_new.setContactStatus(robot);
    robot.generateFeasibleConfiguration(s_new.q);
    s_new.v = Eigen::VectorXd::Random(robot.dimv());
    s_new.a = Eigen::VectorXd::Random(robot.dimv());
    s_new.f = Eigen::VectorXd::Random(robot.max_dimf());
    s_new.mu = Eigen::VectorXd::Random(robot.dim_passive()+robot.max_dimf());
    s_new.lmd = Eigen::VectorXd::Random(robot.dimv());
    s_new.gmm = Eigen::VectorXd::Random(robot.dimv());
    d = SplitDirection(robot);
    d.dq() = Eigen::VectorXd::Random(robot.dimv());
    d.dv() = Eigen::VectorXd::Random(robot.dimv());
    d.da() = Eigen::VectorXd::Random(robot.dimv());
    d.df() = Eigen::VectorXd::Random(robot.dimf());
    d.du = Eigen::VectorXd::Random(robot.dimv());
    d_prev = SplitDirection(robot);
    d_prev.dq() = Eigen::VectorXd::Random(robot.dimv());
    d_prev.dv() = Eigen::VectorXd::Random(robot.dimv());
    d_prev.da() = Eigen::VectorXd::Random(robot.dimv());
    d_prev.df() = Eigen::VectorXd::Random(robot.dimf());
    d_prev.du = Eigen::VectorXd::Random(robot.dimv());
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_prev = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_prev);
    v_prev = Eigen::VectorXd::Random(robot.dimv());
    dq_prev = Eigen::VectorXd::Random(robot.dimv());
    dv_prev = Eigen::VectorXd::Random(robot.dimv());
    std::shared_ptr<JointSpaceCost> joint_cost = std::make_shared<JointSpaceCost>(robot);
    std::shared_ptr<ContactCost> contact_cost = std::make_shared<ContactCost>(robot);
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
  SplitSolution s, s_next, s_tmp, s_old, s_new;
  SplitDirection d, d_prev;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  StateEquation state_equation;
  RobotDynamics robot_dynamics;
  Eigen::VectorXd q_prev, v_prev, dq_prev, dv_prev;
};


TEST_F(FixedBaseTerminalParNMPCTest, isFeasible) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  EXPECT_EQ(parnmpc.isFeasible(robot, s), 
            constraints->isFeasible(robot, constraints_data, s));
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.Fq() = q_prev - s.q + dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  kkt_residual.lq() = - s.lmd;
  kkt_residual.lv() = - s.gmm + dtau * s.lmd;
  kkt_residual.la() = dtau * s.gmm;
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()
                         + kkt_residual.Fv().squaredNorm()
                         + kkt_residual.lq().squaredNorm()
                         + kkt_residual.lv().squaredNorm()
                         + kkt_residual.la().squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = parnmpc.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNormStateEquationAndInverseDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.Fq() = q_prev - s.q + dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  kkt_residual.lq() = - s.lmd;
  kkt_residual.lv() = - s.gmm + dtau * s.lmd;
  kkt_residual.la() = dtau * s.gmm;
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
                         + dtau*dtau*kkt_residual.u_res.squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
  auto pair = parnmpc.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNormStateEquationAndRobotDynamics) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({true}));
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.Fq() = q_prev - s.q + dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  kkt_residual.lq() = - s.lmd;
  kkt_residual.lv() = - s.gmm + dtau * s.lmd;
  kkt_residual.la() = dtau * s.gmm;
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
  auto pair = parnmpc.costAndViolation(robot, t, dtau, s);
  EXPECT_DOUBLE_EQ(pair.first, 0);
  const double violation_ref = kkt_residual.Fx().lpNorm<1>() 
                                + dtau * kkt_residual.u_res.lpNorm<1>()
                                + kkt_residual.C().lpNorm<1>();
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitParNMPC parnmpc(robot, empty_cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual.lu);
  state_equation.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  state_equation.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, KKTErrorNorm) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, 
                                                               q_prev, v_prev, s);
  kkt_residual.setContactStatus(robot);
  kkt_matrix.setContactStatus(robot);
  s.setContactStatus(robot);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual.lu);
  state_equation.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  robot_dynamics.augmentRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm(dtau);
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, costAndViolation) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, s_next);
  const auto pair = parnmpc.costAndViolationTerminal(robot, t, dtau, s); 
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s) 
          + cost->phi(robot, cost_data, t, s)
          + constraints->costSlackBarrier(constraints_data);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  kkt_residual.Fq() = qdiff + dtau * s.v;
  kkt_residual.Fv() = v_prev - s.v + dtau * s.a;
  robot.setContactForces(s.f);
  robot.RNEA(s.q, s.v, s.a, kkt_residual.u_res);
  kkt_residual.u_res -= s.u;
  robot.updateKinematics(s.q, s.v, s.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, costAndViolationWithStepSize) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, s_next);
  const double step_size = 0.3;
  const auto pair = parnmpc.costAndViolationTerminal(robot, step_size, t, dtau, s_old, 
                                                     d_prev, s, d, s_new); 
  s_new.q = s.q + step_size * d.dq();
  s_new.v = s.v + step_size * d.dv();
  s_new.a = s.a + step_size * d.da();
  s_new.f_active() = s.f_active() + step_size * d.df();
  s_new.u = s.u + step_size * d.du;
  const double cost_ref 
      = cost->l(robot, cost_data, t, dtau, s_new) 
          + cost->phi(robot, cost_data, t, s_new)
          + constraints->costSlackBarrier(constraints_data, step_size);
  EXPECT_DOUBLE_EQ(pair.first, cost_ref);
  kkt_residual.Fq() = s_old.q + step_size * d_prev.dq() - s_new.q + dtau * s_new.v;
  kkt_residual.Fv() = s_old.v + step_size * d_prev.dv() - s_new.v + dtau * s_new.a;
  robot.setContactForces(s_new.f);
  robot.RNEA(s_new.q, s_new.v, s_new.a, kkt_residual.u_res);
  kkt_residual.u_res -= s_new.u;
  robot.updateKinematics(s_new.q, s_new.v, s_new.a);
  robot.computeBaumgarteResidual(dtau, kkt_residual.C());
  const double violation_ref 
      = kkt_residual.Fq().lpNorm<1>() + kkt_residual.Fv().lpNorm<1>() 
          + dtau * kkt_residual.u_res.lpNorm<1>() 
          + kkt_residual.C().head(robot.dimf()).lpNorm<1>()
          + constraints->residualL1Nrom(robot, constraints_data, dtau, s_new);
  EXPECT_DOUBLE_EQ(pair.second, violation_ref);
}


TEST_F(FixedBaseTerminalParNMPCTest, coarseUpdate) {
  s.setContactStatus(robot);
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  Eigen::MatrixXd aux_mat_seed = Eigen::MatrixXd::Random(2*robot.dimv(), 
                                                         2*robot.dimv());
  const Eigen::MatrixXd aux_mat_next = aux_mat_seed.transpose() * aux_mat_seed;
  SplitSolution s_new_coarse(robot);
  s_new_coarse.setContactStatus(robot);
  parnmpc.coarseUpdateTerminal(robot, t, dtau, q_prev, v_prev, s, d, s_new_coarse);
  // coarse update ref
  kkt_residual.lu.setZero();
  kkt_matrix.Quu.setZero();
  cost->lu(robot, cost_data, t, dtau, s.u, kkt_residual.lu);
  constraints->augmentDualResidual(robot, constraints_data, dtau, 
                                    kkt_residual.lu);
  cost->luu(robot, cost_data, t, dtau, s.u, kkt_matrix.Quu);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s.u, 
                                     kkt_matrix.Quu, kkt_residual.lu);
  robot_dynamics.condenseRobotDynamics(robot, dtau, s, kkt_matrix, kkt_residual);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, 
                                     kkt_residual);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, 
                                   kkt_residual);
  state_equation.linearizeBackwardEulerTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  cost->computeStageCostHessian(robot, cost_data, t, dtau, s, kkt_matrix);
  cost->computeTerminalCostHessian(robot, cost_data, t, s, kkt_matrix);
  constraints->condenseSlackAndDual(robot, constraints_data, dtau, s, 
                                     kkt_matrix, kkt_residual);
  kkt_matrix.symmetrize(); 
  const int dimKKT = kkt_residual.dimKKT();
  Eigen::MatrixXd kkt_matrix_inverse = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  kkt_matrix.invert(dtau, kkt_matrix_inverse.topLeftCorner(dimKKT, dimKKT));
  SplitDirection d_ref(robot);
  d_ref.split_direction() = kkt_matrix_inverse.topLeftCorner(dimKKT, dimKKT)
                              * kkt_residual.KKT_residual();
  SplitSolution s_new_coarse_ref(robot);
  s_new_coarse_ref.setContactStatus(robot);
  s_new_coarse_ref.lmd = s.lmd - d_ref.dlmd();
  s_new_coarse_ref.gmm = s.gmm - d_ref.dgmm();
  s_new_coarse_ref.mu_active() = s.mu_active() - d_ref.dmu();
  s_new_coarse_ref.a = s.a - d_ref.da();
  s_new_coarse_ref.f_active() = s.f_active() - d_ref.df();
  robot.integrateConfiguration(s.q, d_ref.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v = s.v - d_ref.dv();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  const int dimx = 2*robot.dimv();
  Eigen::MatrixXd aux_mat = Eigen::MatrixXd::Zero(dimx, dimx);
  parnmpc.getAuxiliaryMatrix(aux_mat);
  Eigen::MatrixXd aux_mat_ref = - kkt_matrix_inverse.topLeftCorner(dimx, dimx);
  EXPECT_TRUE(aux_mat.isApprox(aux_mat_ref));
  parnmpc.backwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
  Eigen::VectorXd x_res = Eigen::VectorXd::Zero(dimx);
  x_res.head(robot.dimv()) = s_new.lmd - s_old.lmd;
  x_res.tail(robot.dimv()) = s_new.gmm - s_old.gmm;
  Eigen::VectorXd dx = kkt_matrix_inverse.topRightCorner(dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= dx.head(robot.dimv());
  s_new_coarse_ref.gmm -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  parnmpc.backwardCorrectionParallel(robot, d, s_new_coarse);
  d_ref.split_direction().tail(dimKKT-dimx)
      = kkt_matrix_inverse.bottomRightCorner(dimKKT-dimx, dimx) * x_res;
  s_new_coarse_ref.mu_active().noalias() -= d_ref.dmu();
  s_new_coarse_ref.a.noalias() -= d_ref.da();
  s_new_coarse_ref.f_active().noalias() -= d_ref.df();
  robot.integrateConfiguration(d_ref.dq(), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v.noalias() -= d_ref.dv();
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCorrectionSerial(robot, s_old, s_new, s_new_coarse);
  robot.subtractConfiguration(s_new.q, s_old.q, x_res.head(robot.dimv()));
  x_res.tail(robot.dimv()) = s_new.v - s_old.v;
  dx = kkt_matrix_inverse.bottomLeftCorner(dimx, dimx) * x_res;
  robot.integrateConfiguration(dx.head(robot.dimv()), -1, s_new_coarse_ref.q);
  s_new_coarse_ref.v -= dx.tail(robot.dimv());
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
  parnmpc.forwardCorrectionParallel(robot, d, s_new_coarse);
  d_ref.split_direction().head(dimKKT-dimx) = kkt_matrix_inverse.topLeftCorner(dimKKT-dimx, dimx) * x_res;
  s_new_coarse_ref.lmd -= d_ref.dlmd();
  s_new_coarse_ref.gmm -= d_ref.dgmm();
  s_new_coarse_ref.mu_active() -= d_ref.dmu();
  s_new_coarse_ref.a -= d_ref.da();
  s_new_coarse_ref.f_active() -= d_ref.df();
  EXPECT_TRUE(s_new_coarse.lmd.isApprox(s_new_coarse_ref.lmd));
  EXPECT_TRUE(s_new_coarse.gmm.isApprox(s_new_coarse_ref.gmm));
  EXPECT_TRUE(s_new_coarse.mu.isApprox(s_new_coarse_ref.mu));
  EXPECT_TRUE(s_new_coarse.a.isApprox(s_new_coarse_ref.a));
  EXPECT_TRUE(s_new_coarse.f.isApprox(s_new_coarse_ref.f));
  EXPECT_TRUE(s_new_coarse.q.isApprox(s_new_coarse_ref.q));
  EXPECT_TRUE(s_new_coarse.v.isApprox(s_new_coarse_ref.v));
}


TEST_F(FixedBaseTerminalParNMPCTest, computePrimalDualDirection) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  parnmpc.computePrimalAndDualDirection(robot, dtau, s, s_new, d);
  EXPECT_TRUE(d.dlmd().isApprox(s_new.lmd-s.lmd));
  EXPECT_TRUE(d.dgmm().isApprox(s_new.gmm-s.gmm));
  EXPECT_TRUE(d.dmu().isApprox(s_new.mu_active()-s.mu_active()));
  EXPECT_TRUE(d.da().isApprox(s_new.a-s.a));
  EXPECT_TRUE(d.df().isApprox(s_new.f_active()-s.f_active()));
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s_new.q, s.q, qdiff);
  EXPECT_TRUE(d.dq().isApprox(qdiff));
  EXPECT_TRUE(d.dv().isApprox(s_new.v-s.v));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}