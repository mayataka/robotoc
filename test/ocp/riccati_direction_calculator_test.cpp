#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/cost/joint_space_cost.hpp"
#include "idocp/cost/contact_force_cost.hpp"
#include "idocp/cost/joint_space_impulse_cost.hpp"
#include "idocp/cost/impulse_force_cost.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/riccati_direction_calculator.hpp"

#include "test_helper.hpp"

namespace idocp {

class RiccatiDirectionCalculatorTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_impulse = 5;
    nproc = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
  }

  virtual void TearDown() {
  }


  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, t, dtau;
};



Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


ContactSequence RiccatiDirectionCalculatorTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dtau);
}


void RiccatiDirectionCalculatorTest::test(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto ocp = OCP(robot, cost, constraints, N, max_num_impulse);
  OCPLinearizer linearizer(N, max_num_impulse, nproc);
  std::vector<Robot> robots(nproc, robot);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  linearizer.initConstraints(ocp, ocp_discretizer, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, ocp_discretizer, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  StateConstraintRiccatiFactorizer constraint_factorizer(robot, N, max_num_impulse, nproc);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  if (ocp_discretizer.existImpulse()) {
    riccati_recursion.forwardStateConstraintFactorization(
        factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
    riccati_recursion.backwardStateConstraintFactorization(
        factorizer, ocp_discretizer, kkt_matrix, constraint_factorization);
  }
  const int num_impulse = ocp_discretizer.numImpulseStages();
  const int num_lift = ocp_discretizer.numLiftStages();
  Direction d = Direction(robot, N, max_num_impulse);
  auto d_ref = d;
  RiccatiDirectionCalculator::computeInitialStateDirection(robots, q, v, s, d);
  robot.subtractConfiguration(q, s[0].q, d_ref[0].dq());
  d_ref[0].dv() = v - s[0].v;
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  if (ocp_discretizer.existImpulse()) {
    constraint_factorizer.computeLagrangeMultiplierDirection(
        ocp_discretizer, factorization, constraint_factorization, d);
    constraint_factorizer.aggregateLagrangeMultiplierDirection(
        constraint_factorization, ocp_discretizer, d, factorization);
  }
  riccati_recursion.forwardRiccatiRecursion(
      factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization, d);
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  d_ref = d;
  auto ocp_ref = ocp;
  RiccatiDirectionCalculator direction_calculator(N, max_num_impulse, nproc);
  direction_calculator.computeNewtonDirectionFromRiccatiFactorization(
      ocp, ocp_discretizer, robots, factorizer, factorization, s, d);
  const double primal_step_size = direction_calculator.maxPrimalStepSize();
  const double dual_step_size = direction_calculator.maxDualStepSize();
  const Eigen::VectorXd dx0 = d_ref[0].dx();
  const bool exist_state_constraint = ocp_discretizer.existImpulse();
  double primal_step_size_ref = 1;
  double dual_step_size_ref = 1;
  auto robot_ref = robot;
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_aux = ocp_discretizer.dtau_aux(impulse_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dtau);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                      exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization.impulse[impulse_index], d_ref[i], exist_state_constraint);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dt, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // impulse
      ImpulseSplitRiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d_ref.impulse[impulse_index], 
          exist_state_constraint);
      ocp_ref.impulse[impulse_index].computeCondensedPrimalDirection(
          robot_ref, s.impulse[impulse_index], d_ref.impulse[impulse_index]);
      if (ocp_ref.impulse[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.impulse[impulse_index].maxPrimalStepSize();
      if (ocp_ref.impulse[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.impulse[impulse_index].maxDualStepSize();
      // aux 
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.aux[impulse_index], d_ref.aux[impulse_index], exist_state_constraint);
      factorizer.aux[impulse_index].computeControlInputDirection(
          factorization[i+1], d_ref.aux[impulse_index], exist_state_constraint);
      ocp_ref.aux[impulse_index].computeCondensedPrimalDirection(
          robot_ref, dt_aux, s.aux[impulse_index], d_ref.aux[impulse_index]);
      if (ocp_ref.aux[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.aux[impulse_index].maxPrimalStepSize();
      if (ocp_ref.aux[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.aux[impulse_index].maxDualStepSize();
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_lift = ocp_discretizer.dtau_lift(lift_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dtau);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                      exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization.lift[lift_index], d_ref[i], exist_state_constraint);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dt, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // lift
      SplitRiccatiFactorizer::computeCostateDirection(
          factorization.lift[lift_index], d_ref.lift[lift_index], exist_state_constraint);
      factorizer.lift[lift_index].computeControlInputDirection(
          factorization[i+1], d_ref.lift[lift_index], exist_state_constraint);
      ocp_ref.lift[lift_index].computeCondensedPrimalDirection(
          robot_ref, dt_lift, s.lift[lift_index], d_ref.lift[lift_index]);
      if (ocp_ref.lift[lift_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.lift[lift_index].maxPrimalStepSize();
      if (ocp_ref.lift[lift_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.lift[lift_index].maxDualStepSize();
    }
    else {
      const double dt = dtau;
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], 
                                                      exist_state_constraint);
      factorizer[i].computeControlInputDirection(
          factorization[i+1], d_ref[i], exist_state_constraint);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dt, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
    }
  }
  SplitRiccatiFactorizer::computeCostateDirection(factorization[N], d_ref[N], false);
  if (ocp_ref.terminal.maxPrimalStepSize() < primal_step_size_ref) 
    primal_step_size_ref = ocp_ref.terminal.maxPrimalStepSize();
  if (ocp_ref.terminal.maxDualStepSize() < dual_step_size_ref) 
    dual_step_size_ref = ocp_ref.terminal.maxDualStepSize();
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_DOUBLE_EQ(primal_step_size, primal_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
}


TEST_F(RiccatiDirectionCalculatorTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  test(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(RiccatiDirectionCalculatorTest, floatingBase) {
  Robot robot(floating_base_urdf);
  test(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
