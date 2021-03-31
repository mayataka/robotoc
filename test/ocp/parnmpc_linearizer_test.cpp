#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/backward_correction_solver.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class ParNMPCLinearizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N_ideal = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N_ideal;
  }

  virtual void TearDown() {
  }

  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void testComputeKKTResidual(const Robot& robot) const;
  void testIntegrateSolution(const Robot& robot) const;

  int N_ideal, max_num_impulse, nthreads;
  double T, t, dt;
};



Solution ParNMPCLinearizerTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N_ideal, max_num_impulse);
}


Solution ParNMPCLinearizerTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N_ideal, max_num_impulse, t, true);
}


ContactSequence ParNMPCLinearizerTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N_ideal, max_num_impulse, t, 3*dt);
}


void ParNMPCLinearizerTest::testComputeKKTResidual(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  ParNMPCLinearizer linearizer(N_ideal, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N_ideal, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N_ideal, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot> robots(nthreads, robot);
  auto parnmpc = ParNMPC(robot, cost, constraints, T, N_ideal, max_num_impulse);
  parnmpc.discretize(contact_sequence, t);
  auto parnmpc_ref = parnmpc;
  linearizer.initConstraints(parnmpc, robots, contact_sequence, s);
  linearizer.computeKKTResidual(parnmpc, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  const double kkt_error = linearizer.KKTError(parnmpc, kkt_residual);
  auto robot_ref = robot;
  const int N = parnmpc_ref.discrete().N();
  double kkt_error_ref = 0;
  for (int i=0; i<N; ++i) {
    Eigen::VectorXd q_prev, v_prev;
    if (parnmpc_ref.discrete().isTimeStageAfterImpulse(i)) {
      q_prev = s.impulse[parnmpc_ref.discrete().impulseIndexBeforeTimeStage(i)].q;
      v_prev = s.impulse[parnmpc_ref.discrete().impulseIndexBeforeTimeStage(i)].v;
    }
    else if (parnmpc_ref.discrete().isTimeStageAfterLift(i)) {
      q_prev = s.lift[parnmpc_ref.discrete().liftIndexBeforeTimeStage(i)].q;
      v_prev = s.lift[parnmpc_ref.discrete().liftIndexBeforeTimeStage(i)].v;
    }
    else if (i == 0) {
      q_prev = q;
      v_prev = v;
    }
    else {
      q_prev = s[i-1].q;
      v_prev = s[i-1].v;
    }
    if (i == N-1) {
      const int contact_phase = parnmpc_ref.discrete().contactPhase(i);
      const double dti = parnmpc_ref.discrete().dt(i);
      parnmpc_ref.terminal.initConstraints(robot_ref, i+1, s[i]);
      parnmpc_ref.terminal.computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          parnmpc_ref.discrete().t(i), dti, q_prev, v_prev,
          s[i], kkt_matrix_ref[i], kkt_residual_ref[i]);
      kkt_error_ref += parnmpc_ref.terminal.squaredNormKKTResidual(kkt_residual_ref[i], dt);
    }
    else if (parnmpc_ref.discrete().isTimeStageBeforeImpulse(i)) {
      const int contact_phase = parnmpc_ref.discrete().contactPhase(i);
      const int impulse_index = parnmpc_ref.discrete().impulseIndexAfterTimeStage(i);
      const double ti = parnmpc_ref.discrete().t(i);
      const double t_impulse = parnmpc_ref.discrete().t_impulse(impulse_index);
      const double dti = parnmpc_ref.discrete().dt(i);
      const double dt_aux = parnmpc_ref.discrete().dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      parnmpc_ref[i].initConstraints(robot_ref, i+1, s[i]);
      parnmpc_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, 
          q_prev, v_prev, s[i], s.aux[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      parnmpc_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
      parnmpc_ref.aux[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), t_impulse, dt_aux, 
          s[i].q, s[i].v, s.aux[impulse_index], s.impulse[impulse_index], 
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          contact_sequence.impulseStatus(impulse_index));
      parnmpc_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
      parnmpc_ref.impulse[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, 
          s.aux[impulse_index].q, s.aux[impulse_index].v, 
          s.impulse[impulse_index], s[i+1], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      kkt_error_ref += parnmpc_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
      kkt_error_ref += parnmpc_ref.aux[impulse_index].squaredNormKKTResidual(kkt_residual_ref.aux[impulse_index], dt_aux);
      kkt_error_ref += parnmpc_ref.impulse[impulse_index].squaredNormKKTResidual(kkt_residual_ref.impulse[impulse_index]);
    }
    else if (parnmpc_ref.discrete().isTimeStageBeforeLift(i)) {
      const int contact_phase = parnmpc_ref.discrete().contactPhase(i);
      const int lift_index = parnmpc_ref.discrete().liftIndexAfterTimeStage(i);
      const double ti = parnmpc_ref.discrete().t(i);
      const double t_lift = parnmpc_ref.discrete().t_lift(lift_index);
      const double dti = parnmpc_ref.discrete().dt(i);
      const double dt_lift = parnmpc_ref.discrete().dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      parnmpc_ref[i].initConstraints(robot_ref, i+1, s[i]);
      parnmpc_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, 
          q_prev, v_prev, s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      parnmpc_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
      parnmpc_ref.lift[lift_index].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), t_lift, dt_lift, 
          s[i].q, s[i].v, s.lift[lift_index], s[i+1], 
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
      kkt_error_ref += parnmpc_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
      kkt_error_ref += parnmpc_ref.lift[lift_index].squaredNormKKTResidual(kkt_residual_ref.lift[lift_index], dt_lift);
    }
    else {
      const int contact_phase = parnmpc_ref.discrete().contactPhase(i);
      const double dti = parnmpc_ref.discrete().dt(i);
      parnmpc_ref[i].initConstraints(robot_ref, i+1, s[i]);
      parnmpc_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          parnmpc_ref.discrete().t(i), dti, q_prev, v_prev,
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
      kkt_error_ref += parnmpc_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
    }
  }
  // Discrete events before s[0]. (between s[0] and q, v).
  if (parnmpc_ref.discrete().isTimeStageAfterImpulse(0)) {
    const int contact_phase = parnmpc_ref.discrete().contactPhase(-1);
    const int impulse_index = parnmpc_ref.discrete().impulseIndexBeforeTimeStage(0);
    ASSERT_EQ(contact_phase, 0);
    ASSERT_EQ(impulse_index, 0);
    const double t_impulse = parnmpc_ref.discrete().t_impulse(impulse_index);
    const double dt_aux = parnmpc_ref.discrete().dt_aux(impulse_index);
    parnmpc_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
    parnmpc_ref.aux[impulse_index].computeKKTResidual(
        robot_ref, contact_sequence.contactStatus(contact_phase), 
        t_impulse, dt_aux, q, v, s.aux[impulse_index], s.impulse[impulse_index], 
        kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index]);
    parnmpc_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
    parnmpc_ref.impulse[impulse_index].computeKKTResidual(
        robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, 
        s.aux[impulse_index].q, s.aux[impulse_index].v, 
        s.impulse[impulse_index], s[0], 
        kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
    kkt_error_ref += parnmpc_ref.aux[impulse_index].squaredNormKKTResidual(kkt_residual_ref.aux[impulse_index], dt_aux);
    kkt_error_ref += parnmpc_ref.impulse[impulse_index].squaredNormKKTResidual(kkt_residual_ref.impulse[impulse_index]);
  }
  else if (parnmpc_ref.discrete().isTimeStageAfterLift(0)) {
    const int contact_phase = parnmpc_ref.discrete().contactPhase(-1);
    const int lift_index = parnmpc_ref.discrete().liftIndexBeforeTimeStage(0);
    ASSERT_EQ(contact_phase, 0);
    ASSERT_EQ(lift_index, 0);
    const double t_lift = parnmpc_ref.discrete().t_lift(lift_index);
    const double dt_lift = parnmpc_ref.discrete().dt_lift(lift_index);
    parnmpc_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
    parnmpc_ref.lift[lift_index].computeKKTResidual(
        robot_ref, contact_sequence.contactStatus(contact_phase), t_lift, dt_lift, 
        q, v, s.lift[lift_index], s[0], 
        kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
    kkt_error_ref += parnmpc_ref.lift[lift_index].squaredNormKKTResidual(kkt_residual_ref.lift[lift_index], dt_lift);
  }
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual_ref));
  EXPECT_DOUBLE_EQ(kkt_error, std::sqrt(kkt_error_ref));
}


void ParNMPCLinearizerTest::testIntegrateSolution(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N_ideal, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N_ideal, max_num_impulse);
  auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::vector<Robot> robots(nthreads, robot);
  auto parnmpc = ParNMPC(robot, cost, constraints, T, N_ideal, max_num_impulse);
  parnmpc.discretize(contact_sequence, t);
  ParNMPCLinearizer linearizer(N_ideal, max_num_impulse, nthreads);
  linearizer.initConstraints(parnmpc, robots, contact_sequence, s);
  Direction d(robot, N_ideal, max_num_impulse);
  BackwardCorrection corr(robot, N_ideal, max_num_impulse);
  BackwardCorrectionSolver corr_solver(robot, N_ideal, max_num_impulse, nthreads);
  corr_solver.initAuxMat(parnmpc, robots, s, kkt_matrix);
  corr_solver.coarseUpdate(parnmpc, corr, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  corr_solver.backwardCorrectionSerial(parnmpc, corr, s);
  corr_solver.backwardCorrectionParallel(parnmpc, corr, robots);
  corr_solver.forwardCorrectionSerial(parnmpc, corr, robots, s);
  corr_solver.forwardCorrectionParallel(parnmpc, corr, robots, kkt_matrix, kkt_residual, s, d);
  const double primal_step_size = corr_solver.maxPrimalStepSize();
  const double dual_step_size = corr_solver.maxDualStepSize();
  ASSERT_TRUE(primal_step_size > 0);
  ASSERT_TRUE(primal_step_size <= 1);
  ASSERT_TRUE(dual_step_size > 0);
  ASSERT_TRUE(dual_step_size <= 1);
  auto parnmpc_ref = parnmpc;
  auto corr_ref = corr;
  auto s_ref = s;
  auto d_ref = d;
  linearizer.integrateSolution(parnmpc, robots, kkt_matrix, kkt_residual, 
                               primal_step_size, dual_step_size, d, s);
  auto robot_ref = robot;
  const int N = parnmpc_ref.discrete().N();
  for (int i=0; i<N; ++i) {
    if (i == N-1) {
      parnmpc_ref.terminal.updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      parnmpc_ref.terminal.updateDual(dual_step_size);
    }
    else if (parnmpc_ref.discrete().isTimeStageBeforeImpulse(i)) {
      const int impulse_index = parnmpc_ref.discrete().impulseIndexAfterTimeStage(i);
      parnmpc_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      parnmpc_ref[i].updateDual(dual_step_size);
      parnmpc_ref.aux[impulse_index].updatePrimal(robot_ref, primal_step_size, 
                                                  d_ref.aux[impulse_index], s_ref.aux[impulse_index]);
      parnmpc_ref.aux[impulse_index].updateDual(dual_step_size);
      parnmpc_ref.impulse[impulse_index].updatePrimal(robot_ref, primal_step_size, 
                                                      d_ref.impulse[impulse_index], s_ref.impulse[impulse_index]);
      parnmpc_ref.impulse[impulse_index].updateDual(dual_step_size);
    }
    else if (parnmpc_ref.discrete().isTimeStageBeforeLift(i)) {
      const int lift_index = parnmpc_ref.discrete().liftIndexAfterTimeStage(i);
      parnmpc_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      parnmpc_ref[i].updateDual(dual_step_size);
      parnmpc_ref.lift[lift_index].updatePrimal(robot_ref, primal_step_size, 
                                                d_ref.lift[lift_index], s_ref.lift[lift_index]);
      parnmpc_ref.lift[lift_index].updateDual(dual_step_size);
    }
    else {
      parnmpc_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      parnmpc_ref[i].updateDual(dual_step_size);
    }
  }
  // Discrete events before s[0]. (between s[0] and q, v).
  if (parnmpc_ref.discrete().isTimeStageAfterImpulse(0)) {
    const int impulse_index = parnmpc_ref.discrete().impulseIndexBeforeTimeStage(0);
    ASSERT_EQ(impulse_index, 0);
    parnmpc_ref.aux[impulse_index].updatePrimal(robot_ref, primal_step_size, 
                                                d_ref.aux[impulse_index], s_ref.aux[impulse_index]);
    parnmpc_ref.aux[impulse_index].updateDual(dual_step_size);
    parnmpc_ref.impulse[impulse_index].updatePrimal(robot_ref, primal_step_size, 
                                                    d_ref.impulse[impulse_index], s_ref.impulse[impulse_index]);
    parnmpc_ref.impulse[impulse_index].updateDual(dual_step_size);
  }
  else if (parnmpc_ref.discrete().isTimeStageAfterLift(0)) {
    const int lift_index = parnmpc_ref.discrete().liftIndexBeforeTimeStage(0);
    ASSERT_EQ(lift_index, 0);
    parnmpc_ref.lift[lift_index].updatePrimal(robot_ref, primal_step_size, 
                                              d_ref.lift[lift_index], s_ref.lift[lift_index]);
    parnmpc_ref.lift[lift_index].updateDual(dual_step_size);
  }
  EXPECT_TRUE(testhelper::IsApprox(s, s_ref));
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
}


TEST_F(ParNMPCLinearizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot();
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
  robot = testhelper::CreateFixedBaseRobot(dt);
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
}


TEST_F(ParNMPCLinearizerTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot();
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
  robot = testhelper::CreateFloatingBaseRobot(dt);
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
