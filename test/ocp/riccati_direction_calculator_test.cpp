#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
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
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
  }

  virtual void TearDown() {
  }

  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution RiccatiDirectionCalculatorTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


ContactSequence RiccatiDirectionCalculatorTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


void RiccatiDirectionCalculatorTest::test(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  StateConstraintJacobian jac(robot, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  std::vector<Robot> robots(nthreads, robot);
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual, jac);
  RiccatiRecursion riccati_recursion(robot, N, nthreads);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(ocp.discrete(), kkt_matrix, kkt_residual, jac, factorization);
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  Direction d = Direction(robot, N, max_num_impulse);
  auto d_ref = d;
  RiccatiDirectionCalculator::computeInitialStateDirection(robots, q, v, kkt_matrix, s, d);
  if (robot.hasFloatingBase()) {
    Eigen::VectorXd dq0(Eigen::VectorXd::Zero(robot.dimv()));
    robot.subtractConfiguration(q, s[0].q, dq0);
    d_ref[0].dq() = dq0;
    d_ref[0].dq().head(6) = - kkt_matrix[0].Fqq_prev_inv * dq0.head(6);
    d_ref[0].dv() = v - s[0].v;
  }
  else {
    d_ref[0].dq() = q - s[0].q;
    d_ref[0].dv() = v - s[0].v;
  }
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  riccati_recursion.forwardRiccatiRecursion(ocp.discrete(), kkt_matrix, kkt_residual, d);
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  d_ref = d;
  auto ocp_ref = ocp;
  RiccatiDirectionCalculator direction_calculator(N, max_num_impulse, nthreads);
  direction_calculator.computeNewtonDirection(ocp, robots, factorization, s, d);
  const double primal_step_size = direction_calculator.maxPrimalStepSize();
  const double dual_step_size = direction_calculator.maxDualStepSize();
  const Eigen::VectorXd dx0 = d_ref[0].dx();
  double primal_step_size_ref = 1;
  double dual_step_size_ref = 1;
  auto robot_ref = robot;
  for (int i=0; i<N; ++i) {
    if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
      const double dti = ocp.discrete().dt(i);
      const double dt_aux = ocp.discrete().dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i]);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dti, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // impulse
      ImpulseSplitRiccatiFactorizer::computeCostateDirection(
          factorization.impulse[impulse_index], d_ref.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].computeCondensedPrimalDirection(
          robot_ref, s.impulse[impulse_index], d_ref.impulse[impulse_index]);
      if (ocp_ref.impulse[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.impulse[impulse_index].maxPrimalStepSize();
      if (ocp_ref.impulse[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.impulse[impulse_index].maxDualStepSize();
      // aux 
      SplitRiccatiFactorizer::computeCostateDirection(factorization.aux[impulse_index], d_ref.aux[impulse_index]);
      ocp_ref.aux[impulse_index].computeCondensedPrimalDirection(
          robot_ref, dt_aux, s.aux[impulse_index], d_ref.aux[impulse_index]);
      if (ocp_ref.aux[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.aux[impulse_index].maxPrimalStepSize();
      if (ocp_ref.aux[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.aux[impulse_index].maxDualStepSize();
    }
    else if (ocp.discrete().isTimeStageBeforeLift(i)) {
      const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
      const double dti = ocp.discrete().dt(i);
      const double dt_lift = ocp.discrete().dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i]);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dti, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // lift
      SplitRiccatiFactorizer::computeCostateDirection(factorization.lift[lift_index], d_ref.lift[lift_index]);
      ocp_ref.lift[lift_index].computeCondensedPrimalDirection(
          robot_ref, dt_lift, s.lift[lift_index], d_ref.lift[lift_index]);
      if (ocp_ref.lift[lift_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.lift[lift_index].maxPrimalStepSize();
      if (ocp_ref.lift[lift_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.lift[lift_index].maxDualStepSize();
      const int time_stage_after_lift 
          = ocp.discrete().timeStageAfterLift(lift_index);
      if (ocp.discrete().isTimeStageBeforeImpulse(time_stage_after_lift)) {
        const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(time_stage_after_lift);
        d_ref.lift[lift_index].setImpulseStatus(contact_sequence.impulseStatus(impulse_index));
        SplitRiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.constraint[impulse_index], d_ref.lift[lift_index]);
      }
    }
    else {
      const double dti = dt;
      SplitRiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i]);
      ocp_ref[i].computeCondensedPrimalDirection(robot_ref, dti, s[i], d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
        d_ref[i].setImpulseStatus(contact_sequence.impulseStatus(impulse_index));
        SplitRiccatiFactorizer::computeLagrangeMultiplierDirection(
            factorization.constraint[impulse_index], d_ref[i]);
      }
    }
  }
  SplitRiccatiFactorizer::computeCostateDirection(factorization[ocp_ref.discrete().N()], 
                                                  d_ref[ocp_ref.discrete().N()]);
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
