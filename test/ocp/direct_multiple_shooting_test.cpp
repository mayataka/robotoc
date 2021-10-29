#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class DirectMultipleShootingTest : public ::testing::Test {
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

  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence) const;
  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;

  void test_computeKKTResidual(const Robot& robot) const;
  void test_computeKKTSystem(const Robot& robot) const;
  void test_integrateSolution(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


Solution DirectMultipleShootingTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution DirectMultipleShootingTest::createSolution(const Robot& robot, 
                                                    const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


std::shared_ptr<ContactSequence> DirectMultipleShootingTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequenceSharedPtr(robot, N, max_num_impulse, t, 3*dt);
}


void DirectMultipleShootingTest::test_computeKKTResidual(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  DirectMultipleShooting dms(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  auto ocp_ref = ocp;
  dms.initConstraints(ocp, robots, contact_sequence, s);
  dms.computeKKTResidual(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  const double kkt_error = dms.KKTError(ocp, kkt_residual);
  const double total_cost = dms.totalCost(ocp);
  auto robot_ref = robot;
  double kkt_error_ref = 0;
  double total_cost_ref = 0;
  for (int i=0; i<ocp_ref.discrete().N(); ++i) {
    Eigen::VectorXd q_prev;
    if (i == 0 ) {
      q_prev = q;
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i-1)) {
      q_prev = s.aux[ocp_ref.discrete().impulseIndexAfterTimeStage(i-1)].q;
    }
    else if (ocp_ref.discrete().isTimeStageBeforeLift(i-1)) {
      q_prev = s.lift[ocp_ref.discrete().liftIndexAfterTimeStage(i-1)].q;
    }
    else {
      q_prev = s[i-1].q;
    }
    if (ocp_ref.discrete().isTimeStageBeforeImpulse(i)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(i);
      const double ti = ocp_ref.discrete().t(i);
      const double t_impulse = ocp_ref.discrete().t_impulse(impulse_index);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_aux = ocp_ref.discrete().dt_aux(impulse_index);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.impulse[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence->impulseStatus(impulse_index), t_impulse, s[i].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      ocp_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
      ocp_ref.aux[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase+1), t_impulse, dt_aux, s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[i+1], 
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index]);
      kkt_error_ref += kkt_residual_ref[i].kkt_error;
      kkt_error_ref += kkt_residual_ref.impulse[impulse_index].kkt_error;
      kkt_error_ref += kkt_residual_ref.aux[impulse_index].kkt_error;
      total_cost_ref += ocp_ref[i].stageCost();
      total_cost_ref += ocp_ref.impulse[impulse_index].stageCost();
      total_cost_ref += ocp_ref.aux[impulse_index].stageCost();
    }
    else if (ocp_ref.discrete().isTimeStageBeforeLift(i)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const int lift_index = ocp_ref.discrete().liftIndexAfterTimeStage(i);
      const double ti = ocp_ref.discrete().t(i);
      const double t_lift = ocp_ref.discrete().t_lift(lift_index);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_lift = ocp_ref.discrete().dt_lift(lift_index);
      const bool sto = ocp_ref.discrete().isSTOEnabledLift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
      ocp_ref.lift[lift_index].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase+1), t_lift, dt_lift, s[i].q, 
          s.lift[lift_index], s[i+1], kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
      kkt_error_ref += kkt_residual_ref[i].kkt_error;
      kkt_error_ref += kkt_residual_ref.lift[lift_index].kkt_error;
      total_cost_ref += ocp_ref[i].stageCost();
      total_cost_ref += ocp_ref.lift[lift_index].stageCost();
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i+1)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_next = ocp_ref.discrete().dt(i+1);
      const int impulse_index  
          = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i],
          contact_sequence->impulseStatus(impulse_index),
          dt_next, kkt_matrix_ref.switching[impulse_index], 
          kkt_residual_ref.switching[impulse_index]);
      kkt_error_ref += kkt_residual_ref[i].kkt_error;
      total_cost_ref += ocp_ref[i].stageCost();
    } 
    else {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence->contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
      kkt_error_ref += kkt_residual_ref[i].kkt_error;
      total_cost_ref += ocp_ref[i].stageCost();
    }
  }
  ocp_ref.terminal.computeKKTResidual(robot_ref, t+T, s[ocp_ref.discrete().N()-1].q, 
                                      s[ocp_ref.discrete().N()], 
                                      kkt_matrix_ref[ocp_ref.discrete().N()], 
                                      kkt_residual_ref[ocp_ref.discrete().N()]);
  kkt_error_ref += kkt_residual_ref[ocp_ref.discrete().N()].kkt_error;
  total_cost_ref += ocp_ref.terminal.terminalCost();
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual_ref));
  EXPECT_DOUBLE_EQ(kkt_error, std::sqrt(kkt_error_ref));
  EXPECT_DOUBLE_EQ(total_cost, total_cost_ref);
}


void DirectMultipleShootingTest::test_computeKKTSystem(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  DirectMultipleShooting dms(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  auto ocp_ref = ocp;
  dms.initConstraints(ocp, robots, contact_sequence, s);
  dms.computeKKTSystem(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  auto robot_ref = robot;
  for (int i=0; i<ocp_ref.discrete().N(); ++i) {
    Eigen::VectorXd q_prev;
    if (i == 0) {
      q_prev = q;
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i-1)) {
      q_prev = s.aux[ocp_ref.discrete().impulseIndexAfterTimeStage(i-1)].q;
    }
    else if (ocp_ref.discrete().isTimeStageBeforeLift(i-1)) {
      q_prev = s.lift[ocp_ref.discrete().liftIndexAfterTimeStage(i-1)].q;
    }
    else {
      q_prev = s[i-1].q;
    }
    if (ocp_ref.discrete().isTimeStageBeforeImpulse(i)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(i);
      const double ti = ocp_ref.discrete().t(i);
      const double t_impulse = ocp_ref.discrete().t_impulse(impulse_index);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_aux = ocp_ref.discrete().dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.impulse[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].computeKKTSystem(
          robot_ref, contact_sequence->impulseStatus(impulse_index), t_impulse, s[i].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      ocp_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
      ocp_ref.aux[impulse_index].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase+1), t_impulse, dt_aux, s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[i+1], 
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index]);
    }
    else if (ocp_ref.discrete().isTimeStageBeforeLift(i)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const int lift_index = ocp_ref.discrete().liftIndexAfterTimeStage(i);
      const double ti = ocp_ref.discrete().t(i);
      const double t_lift = ocp_ref.discrete().t_lift(lift_index);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_lift = ocp_ref.discrete().dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
      ocp_ref.lift[lift_index].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase+1), t_lift, dt_lift, s[i].q, 
          s.lift[lift_index], s[i+1], kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i+1)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_next = ocp_ref.discrete().dt(i+1);
      const int impulse_index  
          = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i],
          contact_sequence->impulseStatus(impulse_index),
          dt_next, kkt_matrix_ref.switching[impulse_index],
          kkt_residual_ref.switching[impulse_index]);
    } 
    else {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTSystem(
          robot_ref, contact_sequence->contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
    }
  }
  ocp_ref.terminal.computeKKTSystem(robot_ref, t+T, s[ocp_ref.discrete().N()-1].q, 
                                    s[ocp_ref.discrete().N()], 
                                    kkt_matrix_ref[ocp_ref.discrete().N()], 
                                    kkt_residual_ref[ocp_ref.discrete().N()]);
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual_ref));
}


void DirectMultipleShootingTest::test_integrateSolution(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  DirectMultipleShooting dms(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  dms.initConstraints(ocp, robots, contact_sequence, s);
  dms.computeKKTSystem(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  Direction d(robot, N, max_num_impulse);
  RiccatiRecursion riccati_solver(robots[0], N, max_num_impulse, nthreads);
  RiccatiFactorization riccati_factorization(robots[0], N, max_num_impulse);
  riccati_solver.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, riccati_factorization);
  DirectMultipleShooting::computeInitialStateDirection(ocp, robots, q, v, s, d);
  SplitDirection d0_ref(robot);
  ocp[0].computeInitialStateDirection(robots[0], q, v, s[0], d0_ref);
  EXPECT_TRUE(d[0].isApprox(d0_ref));
  riccati_solver.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  riccati_solver.computeDirection(ocp, riccati_factorization, s, d);
  const double primal_step_size = riccati_solver.maxPrimalStepSize();
  const double dual_step_size = riccati_solver.maxDualStepSize();
  ASSERT_TRUE(primal_step_size > 0);
  ASSERT_TRUE(primal_step_size <= 1);
  ASSERT_TRUE(dual_step_size > 0);
  ASSERT_TRUE(dual_step_size <= 1);
  auto ocp_ref = ocp;
  auto s_ref = s;
  auto d_ref = d;
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  dms.integrateSolution(ocp, robots, primal_step_size, dual_step_size, d, s);
  auto robot_ref = robot;
  for (int i=0; i<ocp_ref.discrete().N(); ++i) {
    if (ocp_ref.discrete().isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_ref.discrete().impulseIndexAfterTimeStage(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_aux = ocp_ref.discrete().dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      const bool sto = ocp.discrete().isSTOEnabledImpulse(impulse_index);
      ocp_ref[i].expandDual(dti, d_ref.impulse[impulse_index], d_ref[i], sto);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
      ocp_ref.impulse[impulse_index].expandDual(d_ref.aux[impulse_index], d_ref.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].updatePrimal(
          robot_ref, primal_step_size, d_ref.impulse[impulse_index], s_ref.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].updateDual(dual_step_size);
      ocp_ref.aux[impulse_index].expandDual(dt_aux, d_ref[i+1], d_ref.aux[impulse_index], sto);
      ocp_ref.aux[impulse_index].updatePrimal(
          robot_ref, primal_step_size, d_ref.aux[impulse_index], s_ref.aux[impulse_index]);
      ocp_ref.aux[impulse_index].updateDual(dual_step_size);
    }
    else if (ocp_ref.discrete().isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_ref.discrete().liftIndexAfterTimeStage(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_lift = ocp_ref.discrete().dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      const bool sto = ocp.discrete().isSTOEnabledLift(lift_index);
      ocp_ref[i].expandDual(dti, d_ref.lift[lift_index], d_ref[i], sto);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
      ocp_ref.lift[lift_index].expandDual(dt_lift, d_ref[i+1], d_ref.lift[lift_index], sto);
      ocp_ref.lift[lift_index].updatePrimal(robot_ref, primal_step_size, d_ref.lift[lift_index], s_ref.lift[lift_index]);
      ocp_ref.lift[lift_index].updateDual(dual_step_size);
    }
    else {
      const double dti = ocp_ref.discrete().dt(i);
      const bool sto_false = false;
      ocp_ref[i].expandDual(dti, d_ref[i+1], d_ref[i], sto_false);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
    }
  }
  ocp_ref.terminal.expandDual(d_ref[ocp_ref.discrete().N()]);
  ocp_ref.terminal.updatePrimal(robot_ref, primal_step_size, d_ref[ocp_ref.discrete().N()], s_ref[ocp_ref.discrete().N()]);
  ocp_ref.terminal.updateDual(dual_step_size);
  EXPECT_TRUE(testhelper::IsApprox(s, s_ref));
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));


  EXPECT_NO_THROW(
    std::cout << s << std::endl;
    std::cout << d << std::endl;
    std::cout << kkt_matrix << std::endl;
    std::cout << kkt_residual << std::endl;
  );
}


TEST_F(DirectMultipleShootingTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot();
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_integrateSolution(robot);
  robot = testhelper::CreateFixedBaseRobot(dt);
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_integrateSolution(robot);
}


TEST_F(DirectMultipleShootingTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot();
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_integrateSolution(robot);
  robot = testhelper::CreateFloatingBaseRobot(dt);
  test_computeKKTResidual(robot);
  test_computeKKTSystem(robot);
  test_integrateSolution(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
