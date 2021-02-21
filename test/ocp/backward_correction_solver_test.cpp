#include <string>
#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"
#include "idocp/ocp/backward_correction_solver.hpp"

#include "test_helper.hpp"

namespace idocp {

class BackwardCorrectionSolverTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
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

  void testCoarseUpdate(const Robot& robot) const;
  void testBackwardCorrection(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N_ideal, max_num_impulse, nthreads;
  double T, t, dt;
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
};



Solution BackwardCorrectionSolverTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N_ideal, max_num_impulse);
}


Solution BackwardCorrectionSolverTest::createSolution(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N_ideal, max_num_impulse, t, true);
}


ContactSequence BackwardCorrectionSolverTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N_ideal, max_num_impulse, t, 3*dt);
}


void BackwardCorrectionSolverTest::testCoarseUpdate(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
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
  ParNMPCLinearizer linearizer(N_ideal, max_num_impulse, nthreads);
  linearizer.initConstraints(parnmpc, robots, contact_sequence, s);
  auto parnmpc_ref = parnmpc;
  BackwardCorrection corr(robot, N_ideal, max_num_impulse);
  auto corr_ref = corr;
  BackwardCorrectionSolver corr_solver(robot, N_ideal, max_num_impulse, nthreads);
  // corr_solver.initAuxMat(parnmpc, robots, s, kkt_matrix);
  corr_solver.coarseUpdate(parnmpc, corr, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  auto robot_ref = robot;
  auto s_new_ref = s;
  const int N = parnmpc_ref.discrete().N();
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  parnmpc_ref.terminal.computeTerminalCostHessian(robot_ref, parnmpc_ref.discrete().t(N-1),
                                                  s[N-1], kkt_matrix_ref[N-1]);
  // const Eigen::MatrixXd aux_mat = kkt_matrix_ref[N-1].Qxx();
  const Eigen::MatrixXd aux_mat = Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv());
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
      parnmpc_ref.terminal.linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          parnmpc_ref.discrete().t(i), dti, q_prev, v_prev,
          s[i], kkt_matrix_ref[i], kkt_residual_ref[i]);
      corr_ref[i].coarseUpdate(robot_ref, dti, kkt_matrix_ref[i], 
                               kkt_residual_ref[i], s[i], s_new_ref[i]);
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
      parnmpc_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, 
          q_prev, v_prev, s[i], s.aux[impulse_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      corr_ref[i].coarseUpdate(robot_ref, dti, aux_mat, kkt_matrix_ref[i], 
                               kkt_residual_ref[i], s[i], s_new_ref[i]);
      parnmpc_ref.aux[impulse_index].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), t_impulse, dt_aux, 
          s[i].q, s[i].v, s.aux[impulse_index], s.impulse[impulse_index], 
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index],
          contact_sequence.impulseStatus(impulse_index));
      corr_ref.aux[impulse_index].coarseUpdate(
          robot_ref, dt_aux, aux_mat, kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], s.aux[impulse_index], 
          s_new_ref.aux[impulse_index]);
      parnmpc_ref.impulse[impulse_index].linearizeOCP(
          robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, 
          s.aux[impulse_index].q, s.aux[impulse_index].v, 
          s.impulse[impulse_index], s[i+1], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
      corr_ref.impulse[impulse_index].coarseUpdate(
          robot_ref, aux_mat, kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], s.impulse[impulse_index], 
          s_new_ref.impulse[impulse_index]);
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
      parnmpc_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), ti, dti, 
          q_prev, v_prev, s[i], s.lift[lift_index], kkt_matrix_ref[i], kkt_residual_ref[i]);
      corr_ref[i].coarseUpdate(robot_ref, dti, aux_mat, kkt_matrix_ref[i], 
                               kkt_residual_ref[i], s[i], s_new_ref[i]);
      parnmpc_ref.lift[lift_index].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), t_lift, dt_lift, 
          s[i].q, s[i].v, s.lift[lift_index], s[i+1], 
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
      corr_ref.lift[lift_index].coarseUpdate(
          robot_ref, dt_lift, aux_mat, kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], s.lift[lift_index], 
          s_new_ref.lift[lift_index]);
    }
    else {
      const int contact_phase = parnmpc_ref.discrete().contactPhase(i);
      const double dti = parnmpc_ref.discrete().dt(i);
      parnmpc_ref[i].linearizeOCP(
          robot_ref, contact_sequence.contactStatus(contact_phase), 
          parnmpc_ref.discrete().t(i), dti, q_prev, v_prev,
          s[i], s[i+1], kkt_matrix_ref[i], kkt_residual_ref[i]);
      corr_ref[i].coarseUpdate(robot_ref, dti, aux_mat, kkt_matrix_ref[i], 
                               kkt_residual_ref[i], s[i], s_new_ref[i]);
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
    parnmpc_ref.aux[impulse_index].linearizeOCP(
        robot_ref, contact_sequence.contactStatus(contact_phase), 
        t_impulse, dt_aux, q, v, s.aux[impulse_index], s.impulse[impulse_index], 
        kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index]);
    corr_ref.aux[impulse_index].coarseUpdate(
        robot_ref, dt_aux, aux_mat, kkt_matrix_ref.aux[impulse_index], 
        kkt_residual_ref.aux[impulse_index], s.aux[impulse_index], 
        s_new_ref.aux[impulse_index]);
    parnmpc_ref.impulse[impulse_index].linearizeOCP(
        robot_ref, contact_sequence.impulseStatus(impulse_index), t_impulse, 
        s.aux[impulse_index].q, s.aux[impulse_index].v, 
        s.impulse[impulse_index], s[0], 
        kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index]);
    corr_ref.impulse[impulse_index].coarseUpdate(
        robot_ref, aux_mat, kkt_matrix_ref.impulse[impulse_index], 
        kkt_residual_ref.impulse[impulse_index], s.impulse[impulse_index], 
        s_new_ref.impulse[impulse_index]);
  }
  else if (parnmpc_ref.discrete().isTimeStageAfterLift(0)) {
    const int contact_phase = parnmpc_ref.discrete().contactPhase(-1);
    const int lift_index = parnmpc_ref.discrete().liftIndexBeforeTimeStage(0);
    ASSERT_EQ(contact_phase, 0);
    ASSERT_EQ(lift_index, 0);
    const double t_lift = parnmpc_ref.discrete().t_lift(lift_index);
    const double dt_lift = parnmpc_ref.discrete().dt_lift(lift_index);
    parnmpc_ref.lift[lift_index].linearizeOCP(
        robot_ref, contact_sequence.contactStatus(contact_phase), t_lift, dt_lift, 
        q, v, s.lift[lift_index], s[0], 
        kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index]);
    corr_ref.lift[lift_index].coarseUpdate(
        robot_ref, dt_lift, aux_mat, kkt_matrix_ref.lift[lift_index], 
        kkt_residual_ref.lift[lift_index], s.lift[lift_index], 
        s_new_ref.lift[lift_index]);
  }
  for (int i=0; i<N; ++i) {
    std::cout << "i = " << i << std::endl;
    if (parnmpc_ref.discrete().isTimeStageAfterImpulse(i)) {
      std::cout << "isTimeStageAfterImpulse!" << std::endl;
    }
    if (parnmpc_ref.discrete().isTimeStageAfterLift(i)) {
      std::cout << "is TimeStageAfterLift!" << std::endl;
    }
    EXPECT_TRUE(kkt_matrix[i].isApprox(kkt_matrix_ref[i]));
  }
  for (int i=0; i<parnmpc_ref.discrete().N_impulse(); ++i) {
    EXPECT_TRUE(kkt_matrix.impulse[i].isApprox(kkt_matrix_ref.impulse[i]));
  }
  for (int i=0; i<parnmpc_ref.discrete().N_impulse(); ++i) {
    EXPECT_TRUE(kkt_matrix.aux[i].isApprox(kkt_matrix_ref.aux[i]));
  }
  for (int i=0; i<parnmpc_ref.discrete().N_lift(); ++i) {
    EXPECT_TRUE(kkt_matrix.lift[i].isApprox(kkt_matrix_ref.lift[i]));
  }
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix_ref));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual_ref));
}


void BackwardCorrectionSolverTest::testBackwardCorrection(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = KKTMatrix(robot, N_ideal, max_num_impulse);
  auto kkt_residual = KKTResidual(robot, N_ideal, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  std::vector<Robot> robots(nthreads, robot);
  auto parnmpc = ParNMPC(robot, cost, constraints, T, N_ideal, max_num_impulse);
  parnmpc.discretize(contact_sequence, t);
  ParNMPCLinearizer linearizer(N_ideal, max_num_impulse, nthreads);
  linearizer.initConstraints(parnmpc, robots, contact_sequence, s);
  auto parnmpc_ref = parnmpc;
  BackwardCorrection corr(robot, N_ideal, max_num_impulse);
  BackwardCorrectionSolver corr_solver(robot, N_ideal, max_num_impulse, nthreads);
  corr_solver.initAuxMat(parnmpc, robots, s, kkt_matrix);
  corr_solver.coarseUpdate(parnmpc, corr, robots, contact_sequence, q, v, s, kkt_matrix, kkt_residual);
  auto corr_ref = corr;
  auto robot_ref = robot;
  auto s_new_ref = s;
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  auto kkt_matrix_ref = kkt_matrix;
  auto kkt_residual_ref = kkt_residual;

  Direction d(robot, N_ideal, max_num_impulse);

  corr_solver.backwardCorrectionSerial(parnmpc, corr, s);
  corr_solver.backwardCorrectionParallel(parnmpc, corr, robots);
  corr_solver.forwardCorrectionSerial(parnmpc, corr, robots, s);
  corr_solver.forwardCorrectionParallel(parnmpc, corr, robots, kkt_matrix, kkt_residual, s, d);
}


TEST_F(BackwardCorrectionSolverTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testCoarseUpdate(robot);
  testBackwardCorrection(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  testCoarseUpdate(robot);
  testBackwardCorrection(robot);
}


TEST_F(BackwardCorrectionSolverTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testCoarseUpdate(robot);
  testBackwardCorrection(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  testCoarseUpdate(robot);
  testBackwardCorrection(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
