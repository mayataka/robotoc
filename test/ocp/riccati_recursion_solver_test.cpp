#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/riccati_recursion_solver.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"

#include "test_helper.hpp"

namespace idocp {

class RiccatiRecursionSolverTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf, {18});
    floating_base_robot = Robot(floating_base_urdf, {14, 24, 34, 44});
    N = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;
  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  KKTMatrix createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const;
  KKTResidual createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const;

  void testComputeInitialStateDirection(const Robot& robot) const;
  void testRiccatiRecursion(const Robot& robot) const;
  void testComputeDirection(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


ContactSequence RiccatiRecursionSolverTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


Solution RiccatiRecursionSolverTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution RiccatiRecursionSolverTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


KKTMatrix RiccatiRecursionSolverTest::createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTMatrix(robot, contact_sequence, N, max_num_impulse);
}


KKTResidual RiccatiRecursionSolverTest::createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTResidual(robot, contact_sequence, N, max_num_impulse);
}


void RiccatiRecursionSolverTest::testComputeInitialStateDirection(const Robot& robot) const {
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  const auto s = createSolution(robot);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::vector<Robot> robots(nthreads, robot);
  Direction d = Direction(robot, N, max_num_impulse);
  auto d_ref = d;
  RiccatiRecursionSolver::computeInitialStateDirection(robots, q, v, kkt_matrix, s, d);
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
}


void RiccatiRecursionSolverTest::testRiccatiRecursion(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  const auto contact_sequence = createContactSequence(robot);
  KKTMatrix kkt_matrix(robot, N, max_num_impulse);
  KKTResidual kkt_residual(robot, N, max_num_impulse);
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  StateConstraintJacobian jac(robot, max_num_impulse);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto s = testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
  linearizer.linearizeOCP(ocp, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual, jac);
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  auto jac_ref = jac;
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  auto factorization_ref = factorization;
  RiccatiRecursionSolver riccati_solver(robot, N, max_num_impulse, nthreads);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  const auto ocp_discretizer = ocp.discrete();
  riccati_solver.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, jac, factorization);
  factorization_ref[ocp_discretizer.N()].Pqq = kkt_matrix_ref[ocp_discretizer.N()].Qqq();
  factorization_ref[ocp_discretizer.N()].Pvv = kkt_matrix_ref[ocp_discretizer.N()].Qvv();
  factorization_ref[ocp_discretizer.N()].sq = - kkt_residual_ref[ocp_discretizer.N()].lq();
  factorization_ref[ocp_discretizer.N()].sv = - kkt_residual_ref[ocp_discretizer.N()].lv();
  for (int i=N-1; i>=0; --i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_aux = ocp_discretizer.dt_aux(impulse_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dt);
      factorizer.aux[impulse_index].backwardRiccatiRecursion(
          factorization_ref[i+1], dt_aux, kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], factorization_ref.aux[impulse_index]);
      factorizer.impulse[impulse_index].backwardRiccatiRecursion(
          factorization_ref.aux[impulse_index], kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], factorization_ref.impulse[impulse_index]);
      factorizer[i].backwardRiccatiRecursion(
          factorization_ref.impulse[impulse_index], dti, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_lift = ocp_discretizer.dt_lift(lift_index);
      ASSERT_TRUE(dti >= 0);
      ASSERT_TRUE(dti <= dt);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dt);
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer.lift[lift_index].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt_lift(lift_index), 
            kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            jac_ref[impulse_index], factorization_ref.lift[lift_index],
            factorization_ref.constraint[impulse_index]);
      }
      else {
        factorizer.lift[lift_index].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt_lift(lift_index), 
            kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
            factorization_ref.lift[lift_index]);
      }
      factorizer[i].backwardRiccatiRecursion(
          factorization_ref.lift[lift_index], dti, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i+1)) {
        const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i+1);
        factorizer[i].backwardRiccatiRecursion(
            factorization_ref[i+1], ocp_discretizer.dt(i), kkt_matrix_ref[i], 
            kkt_residual_ref[i], jac_ref[impulse_index], factorization_ref[i], 
            factorization_ref.constraint[impulse_index]);
      }
      else {
        factorizer[i].backwardRiccatiRecursion(
            factorization_ref[i+1], dt, kkt_matrix_ref[i], 
            kkt_residual_ref[i], factorization_ref[i]);
      }
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.constraint[i].isApprox(factorization_ref.constraint[i]));
  }
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=ocp_discretizer.N(); ++i) {
    d[i].dx().setRandom();
    d[i].du().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].dx().setRandom();
    d.aux[i].dx().setRandom();
    d.aux[i].du().setRandom();
    d.lift[i].dx().setRandom();
    d.lift[i].du().setRandom();
  }
  auto d_ref = d;
  riccati_solver.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  for (int i=0; i<ocp_discretizer.N(); ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_aux = ocp_discretizer.dt_aux(impulse_index);
      factorizer[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], dti, d_ref[i], d_ref.impulse[impulse_index]);
      factorizer.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          d_ref.impulse[impulse_index], d_ref.aux[impulse_index]);
      factorizer.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          dt_aux, d_ref.aux[impulse_index], d_ref[i+1]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndexAfterTimeStage(i);
      const double dti = ocp_discretizer.dt(i);
      const double dt_lift = ocp_discretizer.dt_lift(lift_index);
      factorizer[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          dti, d_ref[i], d_ref.lift[lift_index]);
      factorizer.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          dt_lift, d_ref.lift[lift_index], d_ref[i+1]);
    }
    else {
      factorizer[i].forwardRiccatiRecursion(kkt_matrix_ref[i], 
                                            kkt_residual_ref[i], dt, 
                                            d_ref[i], d_ref[i+1]);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
}


void RiccatiRecursionSolverTest::testComputeDirection(const Robot& robot) const {
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
  RiccatiRecursionSolver riccati_solver(robot, N, max_num_impulse, nthreads);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  riccati_solver.backwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, jac, factorization);
  const int N_impulse = ocp.discrete().N_impulse();
  const int N_lift = ocp.discrete().N_lift();
  Direction d = Direction(robot, N, max_num_impulse);
  auto d_ref = d;
  RiccatiRecursionSolver::computeInitialStateDirection(robots, q, v, kkt_matrix, s, d);
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
  riccati_solver.forwardRiccatiRecursion(ocp, kkt_matrix, kkt_residual, d);
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  d_ref = d;
  auto ocp_ref = ocp;
  riccati_solver.computeDirection(ocp, robots, factorization, s, d);
  const double primal_step_size = riccati_solver.maxPrimalStepSize();
  const double dual_step_size = riccati_solver.maxDualStepSize();
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


TEST_F(RiccatiRecursionSolverTest, fixedBase) {
  testComputeInitialStateDirection(fixed_base_robot);
  testRiccatiRecursion(fixed_base_robot);
  testComputeDirection(fixed_base_robot);
}


TEST_F(RiccatiRecursionSolverTest, floating_base) {
  testComputeInitialStateDirection(floating_base_robot);
  testRiccatiRecursion(floating_base_robot);
  testComputeDirection(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}