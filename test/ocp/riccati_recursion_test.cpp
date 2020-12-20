#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/split_riccati_factorizer.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"

#include "test_helper.hpp"

namespace idocp {

class RiccatiRecursionTest : public ::testing::Test {
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
    nproc = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  KKTMatrix createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const;
  KKTResidual createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void testIsConstraintFactorizationSame(const StateConstraintRiccatiFactorization& lhs, 
                                         const StateConstraintRiccatiFactorization& rhs) const;

  void testBackwardRiccatiRecursion(const Robot& robot) const;
  void testForwardRiccatiRecursionParallel(const Robot& robot) const;
  void testForwardStateConstraintFactorization(const Robot& robot) const;
  void testBackwardStateConstraintFactorization(const Robot& robot) const;
  void testAggregateLagrangeMultiplierDirection(const Robot& robot) const;
  void testForwardRiccatiRecursion(const Robot& robot) const;

  int N, max_num_impulse, nproc;
  double T, t, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


KKTMatrix RiccatiRecursionTest::createKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTMatrix(robot, contact_sequence, N, max_num_impulse);
}


KKTResidual RiccatiRecursionTest::createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const {
  return testhelper::CreateKKTResidual(robot, contact_sequence, N, max_num_impulse);
}


ContactSequence RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dtau);
}


void RiccatiRecursionTest::testIsConstraintFactorizationSame(
    const StateConstraintRiccatiFactorization& lhs, 
    const StateConstraintRiccatiFactorization& rhs) const {
  EXPECT_TRUE(lhs.ENT().isApprox(rhs.ENT()));
  EXPECT_TRUE(lhs.e().isApprox(rhs.e()));
  EXPECT_TRUE(lhs.dxi().isApprox(rhs.dxi()));
  for (int constraint_index=0; constraint_index<max_num_impulse; ++constraint_index) {
    for (int i=0; i<N; ++i) {
      EXPECT_TRUE(lhs.T(constraint_index, i).isApprox(rhs.T(constraint_index, i)));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(lhs.T_impulse(constraint_index, i).isApprox(rhs.T_impulse(constraint_index, i)));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(lhs.T_aux(constraint_index, i).isApprox(rhs.T_aux(constraint_index, i)));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(lhs.T_lift(constraint_index, i).isApprox(rhs.T_lift(constraint_index, i)));
    }
  }
}


void RiccatiRecursionTest::testBackwardRiccatiRecursion(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  EXPECT_TRUE(factorization[N].Pqq.isApprox(kkt_matrix[N].Qqq()));
  EXPECT_TRUE(factorization[N].Pqv.isZero());
  EXPECT_TRUE(factorization[N].Pvq.isZero());
  EXPECT_TRUE(factorization[N].Pvv.isApprox(kkt_matrix[N].Qvv()));
  EXPECT_TRUE(factorization[N].sq.isApprox(-1*kkt_residual[N].lq()));
  EXPECT_TRUE(factorization[N].sv.isApprox(-1*kkt_residual[N].lv()));
  const SplitRiccatiFactorization riccati_factorization_default = SplitRiccatiFactorization(robot);
  for (int i=N-1; i>=0; --i) { EXPECT_TRUE(factorization[i].isApprox(riccati_factorization_default)); }
  for (int i=0; i<max_num_impulse; ++i) { EXPECT_TRUE(factorization.impulse[i].isApprox(riccati_factorization_default)); }
  for (int i=0; i<max_num_impulse; ++i) { EXPECT_TRUE(factorization.aux[i].isApprox(riccati_factorization_default)); }
  for (int i=0; i<max_num_impulse; ++i) { EXPECT_TRUE(factorization.lift[i].isApprox(riccati_factorization_default)); }
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  auto factorizer_ref = factorizer;
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  for (int i=N-1; i>=0; --i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_aux = ocp_discretizer.dtau_aux(impulse_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dtau);
      factorizer_ref.aux[impulse_index].backwardRiccatiRecursion(
          factorization_ref[i+1], dt_aux, kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], factorization_ref.aux[impulse_index]);
      factorizer_ref.impulse[impulse_index].backwardRiccatiRecursion(
          factorization_ref.aux[impulse_index], kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], factorization_ref.impulse[impulse_index]);
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref.impulse[impulse_index], dt, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_lift = ocp_discretizer.dtau_lift(lift_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dtau);
      factorizer_ref.lift[lift_index].backwardRiccatiRecursion(
          factorization_ref[i+1], dt_lift, kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], factorization_ref.lift[lift_index]);
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref.lift[lift_index], dt, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
    else {
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref[i+1], dtau, kkt_matrix_ref[i], 
          kkt_residual_ref[i], factorization_ref[i]);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
}


void RiccatiRecursionTest::testForwardRiccatiRecursionParallel(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  auto constraint_factorization_ref = constraint_factorization;
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, constraint_factorization);
  const bool exist_state_constraint = ocp_discretizer.existStateConstraint();
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndex(i);
      constraint_factorization_ref.Eq(impulse_index) = kkt_matrix_ref.impulse[impulse_index].Pq();
      constraint_factorization_ref.e(impulse_index) = kkt_residual_ref.impulse[impulse_index].P();
      constraint_factorization_ref.T_impulse(impulse_index, impulse_index).topRows(robot.dimv())
          = kkt_matrix_ref.impulse[impulse_index].Pq().transpose();
      constraint_factorization_ref.T_impulse(impulse_index, impulse_index).bottomRows(robot.dimv()).setZero();
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          exist_state_constraint);
      factorizer_ref.aux[impulse_index].forwardRiccatiRecursionParallel(
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          exist_state_constraint);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndex(i);
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardRiccatiRecursionParallel(
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], kkt_residual_ref[i], 
          exist_state_constraint);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  testIsConstraintFactorizationSame(constraint_factorization, constraint_factorization_ref);
}


void RiccatiRecursionTest::testForwardStateConstraintFactorization(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.forwardStateConstraintFactorization(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  const bool exist_state_constraint = ocp_discretizer.existStateConstraint();
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_aux = ocp_discretizer.dtau_aux(impulse_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_aux >= 0);
      ASSERT_TRUE(dt_aux <= dtau);
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i], kkt_matrix_ref[i], kkt_residual_ref[i], 
          dt, factorization_ref.impulse[impulse_index], exist_state_constraint);
      factorizer_ref.impulse[impulse_index].forwardStateConstraintFactorization(
          factorization_ref.impulse[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          factorization_ref.aux[impulse_index], exist_state_constraint);
      factorizer_ref.aux[impulse_index].forwardStateConstraintFactorization(
          factorization_ref.aux[impulse_index],
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          dt_aux, factorization_ref[i+1], exist_state_constraint);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_lift = ocp_discretizer.dtau_lift(lift_index);
      ASSERT_TRUE(dt >= 0);
      ASSERT_TRUE(dt <= dtau);
      ASSERT_TRUE(dt_lift >= 0);
      ASSERT_TRUE(dt_lift <= dtau);
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i], kkt_matrix_ref[i], kkt_residual_ref[i], 
          dt, factorization_ref.lift[lift_index], exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardStateConstraintFactorization(
          factorization_ref.lift[lift_index], kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], dt_lift, factorization_ref[i+1], 
          exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i], kkt_matrix_ref[i], kkt_residual_ref[i], 
          dtau, factorization_ref[i+1], exist_state_constraint);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
}


void RiccatiRecursionTest::testBackwardStateConstraintFactorization(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, constraint_factorization);
  riccati_recursion.forwardStateConstraintFactorization(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  auto constraint_factorization_ref = constraint_factorization;
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.backwardStateConstraintFactorization(factorizer, ocp_discretizer, kkt_matrix, constraint_factorization);
  const int num_constraint = ocp_discretizer.numImpulseStages();
  for (int constraint_index=0; constraint_index<num_constraint; ++constraint_index) {
    ASSERT_FALSE(constraint_factorization.Eq(constraint_index).isZero());
    ASSERT_FALSE(constraint_factorization.e(constraint_index).isZero());
    ASSERT_FALSE(constraint_factorization.T_impulse(constraint_index, constraint_index).isZero());
    ASSERT_TRUE(constraint_factorization.Eq(constraint_index).isApprox(kkt_matrix.impulse[constraint_index].Pq()));
    ASSERT_TRUE(constraint_factorization.e(constraint_index).isApprox(kkt_residual.impulse[constraint_index].P()));
    const int dimv = robot.dimv();
    ASSERT_TRUE((constraint_factorization.T_impulse(constraint_index, constraint_index).topRows(dimv)).isApprox(kkt_matrix.impulse[constraint_index].Pq().transpose()));
    ASSERT_TRUE((constraint_factorization.T_impulse(constraint_index, constraint_index).bottomRows(dimv)).isZero());
  }
  for (int constraint_index=0; constraint_index<num_constraint; ++constraint_index) {
    const int time_stage_before_constraint = ocp_discretizer.timeStageBeforeImpulse(constraint_index);
    const double dtau_constraint = ocp_discretizer.dtau(time_stage_before_constraint);
    ASSERT_TRUE(dtau_constraint >= 0);
    ASSERT_TRUE(dtau_constraint <= dtau);
    factorizer_ref[time_stage_before_constraint].backwardStateConstraintFactorization(
        constraint_factorization_ref.T_impulse(constraint_index, constraint_index), 
        kkt_matrix[time_stage_before_constraint], dtau_constraint, 
        constraint_factorization_ref.T(constraint_index, time_stage_before_constraint));
    for (int i=time_stage_before_constraint-1; i>=0; --i) {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
        const int impulse_index = ocp_discretizer.impulseIndex(i);
        const double dt = ocp_discretizer.dtau(i);
        const double dt_aux = ocp_discretizer.dtau_aux(impulse_index);
        factorizer_ref.aux[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1),
            kkt_matrix_ref.aux[impulse_index], dt_aux, 
            constraint_factorization_ref.T_aux(constraint_index, impulse_index));
        factorizer_ref.impulse[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_aux(constraint_index, impulse_index),
            kkt_matrix_ref.impulse[impulse_index], 
            constraint_factorization_ref.T_impulse(constraint_index, impulse_index));
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_impulse(constraint_index, impulse_index),
            kkt_matrix_ref[i], dt, 
            constraint_factorization_ref.T(constraint_index, i));
      }
      else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
        const int lift_index = ocp_discretizer.liftIndex(i);
        const double dt = ocp_discretizer.dtau(i);
        const double dt_lift = ocp_discretizer.dtau_lift(lift_index);
        factorizer_ref.lift[lift_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1), 
            kkt_matrix_ref.lift[lift_index], dt_lift, 
            constraint_factorization_ref.T_lift(constraint_index, lift_index));
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_lift(constraint_index, lift_index), 
            kkt_matrix_ref[i], dt, constraint_factorization_ref.T(constraint_index, i));
      }
      else {
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1), 
            kkt_matrix_ref[i], dtau, constraint_factorization_ref.T(constraint_index, i));
      }
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
  testIsConstraintFactorizationSame(constraint_factorization, constraint_factorization_ref);
}


void RiccatiRecursionTest::testForwardRiccatiRecursion(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  OCPDiscretizer ocp_discretizer(T, N, max_num_impulse);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(robot, N, max_num_impulse);
  RiccatiRecursion riccati_recursion(robot, N, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(robot, N, max_num_impulse);
  riccati_recursion.backwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, constraint_factorization);
  riccati_recursion.forwardStateConstraintFactorization(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization);
  riccati_recursion.backwardStateConstraintFactorization(factorizer, ocp_discretizer, kkt_matrix, constraint_factorization);
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  Direction d(robot, N, max_num_impulse);
  for (int i=0; i<=N; ++i) {
    d[i].dx().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].dx().setRandom();
    d.aux[i].dx().setRandom();
    d.lift[i].dx().setRandom();
  }
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(factorizer, ocp_discretizer, kkt_matrix, kkt_residual, factorization, d);
  const bool exist_state_constraint = ocp_discretizer.existStateConstraint();
  for (int i=0; i<N; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_index = ocp_discretizer.impulseIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_aux = ocp_discretizer.dtau_aux(impulse_index);
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], factorization_ref.impulse[impulse_index], 
          d_ref[i], dt, d_ref.impulse[impulse_index], exist_state_constraint);
      factorizer_ref.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.impulse[impulse_index], kkt_residual_ref.impulse[impulse_index], 
          d_ref.impulse[impulse_index], d_ref.aux[impulse_index]);
      factorizer_ref.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.aux[impulse_index], kkt_residual_ref.aux[impulse_index], 
          factorization_ref[i+1], d_ref.aux[impulse_index], dt_aux, d_ref[i+1], exist_state_constraint);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_index = ocp_discretizer.liftIndex(i);
      const double dt = ocp_discretizer.dtau(i);
      const double dt_lift = ocp_discretizer.dtau_lift(lift_index);
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], factorization_ref.lift[lift_index], 
          d_ref[i], dt, d_ref.lift[lift_index], exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix_ref.lift[lift_index], kkt_residual_ref.lift[lift_index], 
          factorization_ref[i+1], d_ref.lift[lift_index], dt_lift, d_ref[i+1], exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], kkt_residual_ref[i], factorization_ref[i+1], 
          d_ref[i], dtau, d_ref[i+1], exist_state_constraint);
    }
  }
  EXPECT_TRUE(testhelper::IsApprox(d, d_ref));
  EXPECT_TRUE(testhelper::IsApprox(factorization, factorization_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_matrix, kkt_matrix_ref));
  EXPECT_TRUE(testhelper::IsApprox(kkt_residual, kkt_residual_ref));
  EXPECT_FALSE(testhelper::HasNaN(factorization));
  EXPECT_FALSE(testhelper::HasNaN(kkt_matrix));
  EXPECT_FALSE(testhelper::HasNaN(kkt_residual));
}


TEST_F(RiccatiRecursionTest, fixedBase) {
  testBackwardRiccatiRecursion(fixed_base_robot);
  testForwardRiccatiRecursionParallel(fixed_base_robot);
  testForwardStateConstraintFactorization(fixed_base_robot);
  testBackwardStateConstraintFactorization(fixed_base_robot);
  testForwardRiccatiRecursion(fixed_base_robot);
}


TEST_F(RiccatiRecursionTest, floating_base) {
  testBackwardRiccatiRecursion(floating_base_robot);
  testForwardRiccatiRecursionParallel(floating_base_robot);
  testForwardStateConstraintFactorization(floating_base_robot);
  testBackwardStateConstraintFactorization(floating_base_robot);
  testForwardRiccatiRecursion(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}