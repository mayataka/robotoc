#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_parnmpc.hpp"
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

class FloatingBaseSplitParNMPCTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf = "../urdf/anymal/anymal.urdf";
    std::vector<int> contact_frames = {14, 24, 34, 44};
    const double baum_a = std::abs(Eigen::VectorXd::Random(1)[0]);
    const double baum_b = std::abs(Eigen::VectorXd::Random(1)[0]);
    robot = Robot(urdf, contact_frames, baum_a, baum_b);
    std::random_device rnd;
    std::vector<bool> contact_status;
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
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    q_prev = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_prev);
    v_prev = Eigen::VectorXd::Random(robot.dimv());
    lmd_next = Eigen::VectorXd::Random(robot.dimv());
    gmm_next = Eigen::VectorXd::Random(robot.dimv());
    q_next = Eigen::VectorXd::Random(robot.dimq());
    robot.normalizeConfiguration(q_next);
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
    const Eigen::VectorXd vf_weight = Eigen::VectorXd::Random(robot.dimv());
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
    inverse_dynamics = InverseDynamics(robot);
  }

  virtual void TearDown() {
  }

  double dtau, t;
  std::string urdf;
  Robot robot;
  std::shared_ptr<CostFunction> cost;
  CostFunctionData cost_data;
  std::shared_ptr<Constraints> constraints;
  ConstraintsData constraints_data;
  SplitSolution s, s_tmp, s_old, s_new;
  SplitDirection d;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  StateEquation state_equation;
  InverseDynamics inverse_dynamics;
  Eigen::VectorXd q_prev, v_prev, lmd_next, gmm_next, q_next;
};


TEST_F(FloatingBaseSplitParNMPCTest, isFeasible) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  EXPECT_EQ(parnmpc.isFeasible(robot, s), 
            constraints->isFeasible(robot, constraints_data, s));
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNorm) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormEmptyCostAndEmptyConstraints) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormEmptyCost) {
  auto empty_cost = std::make_shared<CostFunction>();
  SplitParNMPC parnmpc(robot, empty_cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormEmptyConstraints) {
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormOnlyStateEquation) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  robot.setContactStatus(std::vector<bool>({false, false, false, false}));
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  robot.RNEA(s.q, s.v, s.a, s.u);
  s.beta.setZero();
  s.lmd.setZero();
  s.gmm.setZero();
  s.mu.setZero();
  lmd_next.setZero();
  gmm_next.setZero();
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error = parnmpc.squaredKKTErrorNorm(robot, t, dtau, q_prev, 
                                                       v_prev, s, lmd_next, 
                                                       gmm_next, q_next);
  state_equation.linearizeStateEquation(robot, dtau, q_prev, v_prev, s, 
                                        lmd_next, gmm_next, q_next, 
                                        kkt_matrix, kkt_residual);
  double kkt_error_ref = kkt_residual.Fq().squaredNorm()+kkt_residual.Fv().squaredNorm()+dtau*dtau*s.u.head(6).squaredNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormTerminal) {
  SplitParNMPC parnmpc(robot, cost, constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error 
      = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, q_prev, v_prev, s);
  robot.updateKinematics(s.q, s.v, s.a);
  cost->computeStageCostDerivatives(robot, cost_data, t, dtau, s, kkt_residual);
  cost->computeTerminalCostDerivatives(robot, cost_data, t, s, kkt_residual);
  constraints->augmentDualResidual(robot, constraints_data, dtau, kkt_residual);
  state_equation.linearizeStateEquationTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  kkt_error_ref += constraints->squaredKKTErrorNorm(robot, constraints_data, dtau, s);
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


TEST_F(FloatingBaseSplitParNMPCTest, KKTErrorNormTerminalEmptyCostAndEmptyConstraints) {
  auto empty_cost = std::make_shared<CostFunction>();
  auto empty_constraints = std::make_shared<Constraints>();
  SplitParNMPC parnmpc(robot, empty_cost, empty_constraints);
  parnmpc.initConstraints(robot, 2, dtau, s);
  constraints->setSlackAndDual(robot, constraints_data, dtau, s);
  const double kkt_error 
      = parnmpc.squaredKKTErrorNormTerminal(robot, t, dtau, q_prev, v_prev, s);
  robot.updateKinematics(s.q, s.v, s.a);
  state_equation.linearizeStateEquationTerminal(robot, dtau, q_prev, v_prev, s, 
                                                kkt_matrix, kkt_residual);
  equalityconstraints::LinearizeEqualityConstraints(robot, dtau, s, 
                                                    kkt_matrix, kkt_residual);
  inverse_dynamics.linearizeInverseDynamics(robot, dtau, s, kkt_residual);
  double kkt_error_ref = kkt_residual.squaredKKTErrorNorm();
  EXPECT_DOUBLE_EQ(kkt_error, kkt_error_ref);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}