#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/hybrid/hybrid_container.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "kkt_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class RiccatiRecursionTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
  }

  virtual void TearDown() {
  }

  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;
  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence) const;
  KKTMatrix createKKTMatrix(const Robot& robot, 
                            const std::shared_ptr<ContactSequence>& contact_sequence) const;
  KKTResidual createKKTResidual(const Robot& robot, 
                                const std::shared_ptr<ContactSequence>& contact_sequence) const;
  void test_riccatiRecursion(const Robot& robot) const;
  void test_computeDirection(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


std::shared_ptr<ContactSequence > RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequenceSharedPtr(robot, N, max_num_impulse, t, 3*dt);
}


Solution RiccatiRecursionTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution RiccatiRecursionTest::createSolution(const Robot& robot, 
                                              const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


KKTMatrix RiccatiRecursionTest::createKKTMatrix(const Robot& robot, 
                                                const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateKKTMatrix(robot, contact_sequence, N, max_num_impulse);
}


KKTResidual RiccatiRecursionTest::createKKTResidual(const Robot& robot, 
                                                    const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateKKTResidual(robot, contact_sequence, N, max_num_impulse);
}


TEST_P(RiccatiRecursionTest, riccatiRecursion) {
  const auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  DirectMultipleShooting dms(nthreads);
  const auto contact_sequence = createContactSequence(robot);
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  aligned_vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, contact_sequence, cost, constraints, T, N);
  ocp.discretize(contact_sequence, t);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto s = testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
  dms.initConstraints(ocp, robots, contact_sequence, s);
  dms.computeKKTSystem(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  auto factorization_ref = factorization;
  RiccatiRecursion riccati_recursion(ocp, nthreads);
  RiccatiFactorizer factorizer(robot);
  hybrid_container<LQRPolicy> lqr_policy(robot, N, max_num_impulse);
  const auto discretization = ocp.discrete();
  riccati_recursion.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, factorization);

  factorization_ref[discretization.N()].P = kkt_matrix_ref[discretization.N()].Qxx;
  factorization_ref[discretization.N()].s = - kkt_residual_ref[discretization.N()].lx;
  for (int i=N-1; i>=0; --i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = discretization.impulseIndexAfterTimeStage(i);
      factorizer.backwardRiccatiRecursion(
          factorization_ref[i+1], kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], factorization_ref.aux[impulse_index],
          lqr_policy.aux[impulse_index]);
      factorizer.backwardRiccatiRecursion(
          factorization_ref.aux[impulse_index], kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], factorization_ref.impulse[impulse_index]);
      factorizer.backwardRiccatiRecursion(
          factorization_ref.impulse[impulse_index], kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i], lqr_policy[i]);
    }
    else if (discretization.isTimeStageBeforeLift(i)) {
      const int lift_index = discretization.liftIndexAfterTimeStage(i);
      if (discretization.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = discretization.impulseIndexAfterTimeStage(i+1);
        factorizer.backwardRiccatiRecursion(
            factorization_ref[i+1], kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            kkt_matrix_ref.switching[impulse_index], kkt_residual_ref.switching[impulse_index], 
            factorization_ref.lift[lift_index], factorization_ref.switching[impulse_index], lqr_policy.lift[lift_index]);
      }
      else {
        factorizer.backwardRiccatiRecursion(
            factorization_ref[i+1], 
            kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            factorization_ref.lift[lift_index], lqr_policy.lift[lift_index]);
      }
      factorizer.backwardRiccatiRecursion(
          factorization_ref.lift[lift_index], kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i], lqr_policy[i]);
    }
    else {
      if (discretization.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = discretization.impulseIndexAfterTimeStage(i+1);
        factorizer.backwardRiccatiRecursion(
            factorization_ref[i+1], kkt_matrix_ref[i], kkt_residual_ref[i], 
            kkt_matrix_ref.switching[impulse_index], kkt_residual_ref.switching[impulse_index], 
            factorization_ref[i], factorization_ref.switching[impulse_index], lqr_policy[i]);
      }
      else {
        factorizer.backwardRiccatiRecursion(
            factorization_ref[i+1], kkt_matrix_ref[i], kkt_residual_ref[i], 
            factorization_ref[i], lqr_policy[i]);
      }
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=discretization.N(); ++i) {
    d[i].dx.setRandom();
    d[i].du.setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].dx.setRandom();
    d.aux[i].dx.setRandom();
    d.aux[i].du.setRandom();
    d.lift[i].dx.setRandom();
    d.lift[i].du.setRandom();
  }
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  for (int i=0; i<discretization.N(); ++i) {
    if (discretization.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = discretization.impulseIndexAfterTimeStage(i);
      const bool sto = false;
      const bool sto_next = false;
      factorizer.forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          lqr_policy[i], d_ref[i], d_ref.impulse[impulse_index], sto, sto_next);
      factorizer.forwardRiccatiRecursion(
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          d_ref.impulse[impulse_index], d_ref.aux[impulse_index]);
      const bool sto_next_next = false;
      factorizer.forwardRiccatiRecursion(
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          lqr_policy.aux[impulse_index], d_ref.aux[impulse_index], d_ref[i+1],
          sto_next, sto_next_next);
    }
    else if (discretization.isTimeStageBeforeLift(i)) {
      const int lift_index = discretization.liftIndexAfterTimeStage(i);
      const bool sto = false;
      const bool sto_next = false;
      factorizer.forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          lqr_policy[i], d_ref[i], d_ref.lift[lift_index], sto, sto_next);
      const bool sto_next_next = false;
      factorizer.forwardRiccatiRecursion(
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          lqr_policy.lift[lift_index], d_ref.lift[lift_index], d_ref[i+1],
          sto_next, sto_next_next);
    }
    else {
      const bool sto = false;
      const bool sto_next = false;
      factorizer.forwardRiccatiRecursion(kkt_matrix_ref[i], kkt_residual_ref[i],  
                                         lqr_policy[i], d_ref[i], d_ref[i+1],
                                         sto, sto_next);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  Eigen::MatrixXd Kq(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())), 
                  Kv(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv()));
  for (int i=0; i<N; ++i) {
    riccati_recursion.getStateFeedbackGain(i, Kq, Kv);
    EXPECT_TRUE(lqr_policy[i].Kq().isApprox(Kq));
    EXPECT_TRUE(lqr_policy[i].Kv().isApprox(Kv));
  }
}


TEST_P(RiccatiRecursionTest, computeDirection) {
  const auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto ocp = OCP(robot, contact_sequence, cost, constraints, T, N);
  ocp.discretize(contact_sequence, t);
  DirectMultipleShooting dms(nthreads);
  aligned_vector<Robot> robots(nthreads, robot);
  dms.initConstraints(ocp, robots, contact_sequence, s);
  dms.computeKKTSystem(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  RiccatiRecursion riccati_recursion(ocp, nthreads);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, factorization);
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  Direction d = Direction(robot, N, max_num_impulse);
  dms.computeInitialStateDirection(ocp, robots, q, v, s, d);
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  d_ref = d;
  auto ocp_ref = ocp;
  riccati_recursion.computeDirection(ocp, contact_sequence, factorization, d);
  const double primal_step_size = riccati_recursion.maxPrimalStepSize();
  const double dual_step_size = riccati_recursion.maxDualStepSize();
  const Eigen::VectorXd dx0 = d_ref[0].dx;
  double primal_step_size_ref = 1;
  double dual_step_size_ref = 1;
  auto robot_ref = robot;
  for (int i=0; i<N; ++i) {
    if (ocp.discrete().isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp.discrete().impulseIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      const bool sto_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], sto, sto_next);
      ocp_ref[i].expandPrimal(
          contact_sequence->contactStatus(ocp.discrete().contactPhase(i)), d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // impulse
      RiccatiFactorizer::computeCostateDirection(factorization.impulse[impulse_index], 
                                                 d_ref.impulse[impulse_index], sto);
      ocp_ref.impulse[impulse_index].expandPrimal(
          contact_sequence->impulseStatus(impulse_index), 
          d_ref.impulse[impulse_index]);
      if (ocp_ref.impulse[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.impulse[impulse_index].maxPrimalStepSize();
      if (ocp_ref.impulse[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.impulse[impulse_index].maxDualStepSize();
      // aux 
      const bool sto_next_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization.aux[impulse_index], 
                                                 d_ref.aux[impulse_index], 
                                                 sto_next, sto_next_next);
      ocp_ref.aux[impulse_index].expandPrimal(
          contact_sequence->contactStatus(ocp.discrete().contactPhase(i)+1), 
          d_ref.aux[impulse_index]);
      if (ocp_ref.aux[impulse_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.aux[impulse_index].maxPrimalStepSize();
      if (ocp_ref.aux[impulse_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.aux[impulse_index].maxDualStepSize();
    }
    else if (ocp.discrete().isTimeStageBeforeLift(i)) {
      const int lift_index = ocp.discrete().liftIndexAfterTimeStage(i);
      const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
      const bool sto_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], sto, sto_next);
      ocp_ref[i].expandPrimal(
          contact_sequence->contactStatus(ocp.discrete().contactPhase(i)), d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      // lift
      const bool sto_next_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization.lift[lift_index], 
                                                 d_ref.lift[lift_index], sto_next, sto_next_next);
      ocp_ref.lift[lift_index].expandPrimal(
          contact_sequence->contactStatus(
              ocp.discrete().contactPhaseAfterLift(lift_index)), 
          d_ref.lift[lift_index]);
      if (ocp_ref.lift[lift_index].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref.lift[lift_index].maxPrimalStepSize();
      if (ocp_ref.lift[lift_index].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref.lift[lift_index].maxDualStepSize();
      const int time_stage_after_lift 
          = ocp.discrete().timeStageAfterLift(lift_index);
      ASSERT_FALSE(ocp.discrete().isTimeStageBeforeImpulse(time_stage_after_lift));
    }
    else {
      const bool sto = false;
      const bool sto_next = false;
      RiccatiFactorizer::computeCostateDirection(factorization[i], d_ref[i], sto, sto_next);
      ocp_ref[i].expandPrimal(
          contact_sequence->contactStatus(ocp.discrete().contactPhase(i)), d_ref[i]);
      if (ocp_ref[i].maxPrimalStepSize() < primal_step_size_ref) 
        primal_step_size_ref = ocp_ref[i].maxPrimalStepSize();
      if (ocp_ref[i].maxDualStepSize() < dual_step_size_ref) 
        dual_step_size_ref = ocp_ref[i].maxDualStepSize();
      if (ocp.discrete().isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
        const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
        const bool sto_next = false;
        d_ref[i].setImpulseStatus(contact_sequence->impulseStatus(impulse_index));
        RiccatiFactorizer::computeLagrangeMultiplierDirection(factorization.switching[impulse_index], 
                                                              d_ref[i], sto, sto_next);
      }
    }
  }
  RiccatiFactorizer::computeCostateDirection(factorization[ocp_ref.discrete().N()], 
                                             d_ref[ocp_ref.discrete().N()], false, false);
  if (ocp_ref.terminal.maxPrimalStepSize() < primal_step_size_ref) 
    primal_step_size_ref = ocp_ref.terminal.maxPrimalStepSize();
  if (ocp_ref.terminal.maxDualStepSize() < dual_step_size_ref) 
    dual_step_size_ref = ocp_ref.terminal.maxDualStepSize();
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_DOUBLE_EQ(primal_step_size, primal_step_size_ref);
  EXPECT_DOUBLE_EQ(dual_step_size, dual_step_size_ref);
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, RiccatiRecursionTest, 
  ::testing::Values(testhelper::CreateFixedBaseRobot(),
                    testhelper::CreateFixedBaseRobot(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateFloatingBaseRobot(),
                    testhelper::CreateFloatingBaseRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}