#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"
#include "idocp/ocp/riccati_recursion_solver.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class OCPLinearizerTest : public ::testing::Test {
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
                          const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void testLinearizeOCP(const Robot& robot) const;
  void testComputeKKTResidual(const Robot& robot) const;
  void testIntegrateSolution(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


Solution OCPLinearizerTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution OCPLinearizerTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


ContactSequence OCPLinearizerTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


void OCPLinearizerTest::testLinearizeOCP(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  auto ocp_ref = ocp;
  StateConstraintJacobian jac(robot, max_num_impulse);
  StateConstraintJacobian jac_ref = jac;
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual, jac);
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
      ocp_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.impulse[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].linearizeOCP(
          robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, s[i].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      ocp_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
      ocp_ref.aux[impulse_index].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase+1), t_impulse, dt_aux, s.impulse[impulse_index].q, 
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
      ocp_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
      ocp_ref.lift[lift_index].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase+1), t_lift, dt_lift, s[i].q, 
          s.lift[lift_index], s[i+1], kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i+1)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_next = ocp_ref.discrete().dt(i+1);
      const int impulse_index  
          = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i],
          contact_sequence.impulseStatus(impulse_index),
          dt_next, jac_ref[impulse_index]);
    } 
    else {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
    }
  }
  ocp_ref.terminal.linearizeOCP(robot_ref, t+T, s[ocp_ref.discrete().N()-1].q, 
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


void OCPLinearizerTest::testComputeKKTResidual(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  auto ocp_ref = ocp;
  StateConstraintJacobian jac(robot, max_num_impulse);
  StateConstraintJacobian jac_ref = jac;
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.computeKKTResidual(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual, jac);
  const double kkt_error = linearizer.KKTError(ocp, kkt_residual);
  auto robot_ref = robot;
  double kkt_error_ref = 0;
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
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.impulse[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.impulse[impulse_index].initConstraints(robot_ref, s.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, s[i].q, 
          s.impulse[impulse_index], s.aux[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      ocp_ref.aux[impulse_index].initConstraints(robot_ref, 0, s.aux[impulse_index]);
      ocp_ref.aux[impulse_index].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase+1), t_impulse, dt_aux, s.impulse[impulse_index].q, 
          s.aux[impulse_index], s[i+1], 
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index]);
      kkt_error_ref += ocp_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
      kkt_error_ref += ocp_ref.impulse[impulse_index].squaredNormKKTResidual(kkt_residual_ref.impulse[impulse_index]);
      kkt_error_ref += ocp_ref.aux[impulse_index].squaredNormKKTResidual(kkt_residual_ref.aux[impulse_index], dt_aux);
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
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, q_prev, 
          s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      ocp_ref.lift[lift_index].initConstraints(robot_ref, 0, s.lift[lift_index]);
      ocp_ref.lift[lift_index].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase+1), t_lift, dt_lift, s[i].q, 
          s.lift[lift_index], s[i+1], kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
      kkt_error_ref += ocp_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
      kkt_error_ref += ocp_ref.lift[lift_index].squaredNormKKTResidual(kkt_residual_ref.lift[lift_index], dt_lift);
    }
    else if (ocp_ref.discrete().isTimeStageBeforeImpulse(i+1)) {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      const double dt_next = ocp_ref.discrete().dt(i+1);
      const int impulse_index  
          = ocp_ref.discrete().impulseIndexAfterTimeStage(i+1);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i],
          contact_sequence.impulseStatus(impulse_index),
          dt_next, jac_ref[impulse_index]);
      kkt_error_ref += ocp_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
    } 
    else {
      const int contact_phase = ocp_ref.discrete().contactPhase(i);
      const double dti = ocp_ref.discrete().dt(i);
      ocp_ref[i].initConstraints(robot_ref, i, s[i]);
      ocp_ref[i].computeKKTResidual(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          ocp_ref.discrete().t(i), dti, q_prev, 
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
      kkt_error_ref += ocp_ref[i].squaredNormKKTResidual(kkt_residual_ref[i], dti);
    }
  }
  ocp_ref.terminal.computeKKTResidual(robot_ref, t+T, s[ocp_ref.discrete().N()-1].q, 
                                      s[ocp_ref.discrete().N()], 
                                      kkt_matrix_ref[ocp_ref.discrete().N()], 
                                      kkt_residual_ref[ocp_ref.discrete().N()]);
  kkt_error_ref += ocp_ref.terminal.squaredNormKKTResidual(kkt_residual_ref[ocp_ref.discrete().N()]);
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual_ref));
  EXPECT_DOUBLE_EQ(kkt_error, std::sqrt(kkt_error_ref));
}


void OCPLinearizerTest::testIntegrateSolution(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  StateConstraintJacobian jac(robot, max_num_impulse);
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, robots, contact_sequence, 
                          q, v, s, kkt_matrix, kkt_residual, jac);
  Direction d(robot, N, max_num_impulse);
  RiccatiRecursionSolver riccati_solver(robots[0], N, max_num_impulse, nthreads);
  RiccatiFactorization riccati_factorization(robots[0], N, max_num_impulse);
  riccati_solver.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, 
                                          jac, riccati_factorization);
  riccati_solver.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  riccati_solver.computeDirection(ocp, robots, riccati_factorization, s, d);
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
  linearizer.integrateSolution(ocp, robots, kkt_matrix, kkt_residual, 
                               primal_step_size, dual_step_size, d, s);
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
      ocp_ref[i].computeCondensedDualDirection(
          robot_ref, dti, kkt_matrix_ref[i], kkt_residual_ref[i], 
          d_ref.impulse[impulse_index], d_ref[i]);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
      ocp_ref.impulse[impulse_index].computeCondensedDualDirection(
          robot_ref, kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          d_ref.aux[impulse_index], d_ref.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].updatePrimal(
          robot_ref, primal_step_size, d_ref.impulse[impulse_index], s_ref.impulse[impulse_index]);
      ocp_ref.impulse[impulse_index].updateDual(dual_step_size);
      ocp_ref.aux[impulse_index].computeCondensedDualDirection(
          robot_ref, dt_aux, kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index],
          d_ref[i+1], d_ref.aux[impulse_index]);
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
      ocp_ref[i].computeCondensedDualDirection(
          robot_ref, dti, kkt_matrix_ref[i], kkt_residual_ref[i], 
          d_ref.lift[lift_index], d_ref[i]);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
      ocp_ref.lift[lift_index].computeCondensedDualDirection(
          robot_ref, dt_lift, kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          d_ref[i+1], d_ref.lift[lift_index]);
      ocp_ref.lift[lift_index].updatePrimal(robot_ref, primal_step_size, d_ref.lift[lift_index], s_ref.lift[lift_index]);
      ocp_ref.lift[lift_index].updateDual(dual_step_size);
    }
    else {
      const double dti = ocp_ref.discrete().dt(i);
      ocp_ref[i].computeCondensedDualDirection(
          robot_ref, dti, kkt_matrix_ref[i], kkt_residual_ref[i], d_ref[i+1], d_ref[i]);
      ocp_ref[i].updatePrimal(robot_ref, primal_step_size, d_ref[i], s_ref[i]);
      ocp_ref[i].updateDual(dual_step_size);
    }
  }
  ocp_ref.terminal.computeCondensedDualDirection(robot_ref, 
                                                 kkt_matrix_ref[ocp_ref.discrete().N()], 
                                                 kkt_residual_ref[ocp_ref.discrete().N()], 
                                                 d_ref[ocp_ref.discrete().N()]);
  ocp_ref.terminal.updatePrimal(robot_ref, primal_step_size, d_ref[ocp_ref.discrete().N()], s_ref[ocp_ref.discrete().N()]);
  ocp_ref.terminal.updateDual(dual_step_size);
  EXPECT_TRUE(testhelper::IsApprox(s, s_ref));
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
}


TEST_F(OCPLinearizerTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot();
  testLinearizeOCP(robot);
  // testComputeKKTResidual(robot);
  // testIntegrateSolution(robot);
  robot = testhelper::CreateFixedBaseRobot(dt);
  testLinearizeOCP(robot);
  // testComputeKKTResidual(robot);
  // testIntegrateSolution(robot);
}


TEST_F(OCPLinearizerTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot();
  testLinearizeOCP(robot);
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
  robot = testhelper::CreateFloatingBaseRobot(dt);
  testLinearizeOCP(robot);
  testComputeKKTResidual(robot);
  testIntegrateSolution(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
