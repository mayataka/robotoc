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

  template <typename T>
  void testIsSame(const T& rhs, const T& lhs) const {
    for (int i=0; i<=N; ++i) {
      EXPECT_TRUE(rhs[i].isApprox(lhs[i]));
      EXPECT_FALSE(rhs[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.impulse[i].isApprox(lhs.impulse[i]));
      EXPECT_FALSE(rhs.impulse[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.aux[i].isApprox(lhs.aux[i]));
      EXPECT_FALSE(rhs.aux[i].hasNaN());
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.lift[i].isApprox(lhs.lift[i]));
      EXPECT_FALSE(rhs.lift[i].hasNaN());
    }
  }

  template <typename T>
  void testIsApprox(const T& rhs, const T& lhs) const {
    for (int i=0; i<=N; ++i) {
      EXPECT_TRUE(rhs[i].isApprox(lhs[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.impulse[i].isApprox(lhs.impulse[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.aux[i].isApprox(lhs.aux[i]));
    }
    for (int i=0; i<max_num_impulse; ++i) {
      EXPECT_TRUE(rhs.lift[i].isApprox(lhs.lift[i]));
    }
  }

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
  KKTMatrix kkt_matrix = KKTMatrix(N, max_num_impulse, robot);
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  for (int i=0; i<=N; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.has_floating_base()) {
      kkt_matrix[i].Fqq().setIdentity();
      kkt_matrix[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix[i].Fvq().setRandom();
    kkt_matrix[i].Fvv().setRandom();
    kkt_matrix[i].Fvu().setRandom();
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<num_impulse; ++i) {
    kkt_matrix.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx, dimx);
    kkt_matrix.impulse[i].Qxx() = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    if (robot.has_floating_base()) {
      kkt_matrix.impulse[i].Fqq().setIdentity();
      kkt_matrix.impulse[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.impulse[i].Fvq().setRandom();
    kkt_matrix.impulse[i].Fvv().setRandom();
    kkt_matrix.impulse[i].Pq().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.has_floating_base()) {
      kkt_matrix.aux[i].Fqq().setIdentity();
      kkt_matrix.aux[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.aux[i].Fvq().setRandom();
    kkt_matrix.aux[i].Fvv().setRandom();
    kkt_matrix.aux[i].Fvu().setRandom();
  }
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_lift; ++i) {
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Random(dimx+dimu, dimx+dimu);
    const Eigen::MatrixXd Qxxuu = tmp * tmp.transpose() + Eigen::MatrixXd::Identity(dimx+dimu, dimx+dimu);
    kkt_matrix[i].Qxx() = Qxxuu.topLeftCorner(dimx, dimx);
    kkt_matrix[i].Quu() = Qxxuu.bottomRightCorner(dimu, dimu);
    kkt_matrix[i].Qxu() = Qxxuu.topRightCorner(dimx, dimu);
    if (robot.has_floating_base()) {
      kkt_matrix.lift[i].Fqq().setIdentity();
      kkt_matrix.lift[i].Fqq().topLeftCorner(6, 6).setRandom();
    }
    kkt_matrix.lift[i].Fvq().setRandom();
    kkt_matrix.lift[i].Fvv().setRandom();
    kkt_matrix.lift[i].Fvu().setRandom();
  }
  return kkt_matrix;
}


KKTResidual RiccatiRecursionTest::createKKTResidual(const Robot& robot, const ContactSequence& contact_sequence) const {
  KKTResidual kkt_residual = KKTResidual(N, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    kkt_residual[i].lx().setRandom();
    kkt_residual[i].lu().setRandom();
    kkt_residual[i].Fx().setRandom();
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    kkt_residual.impulse[i].lx().setRandom();
    kkt_residual.impulse[i].Fx().setRandom();
    kkt_residual.impulse[i].P().setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    kkt_residual.aux[i].lx().setRandom();
    kkt_residual.aux[i].lu().setRandom();
    kkt_residual.aux[i].Fx().setRandom();
  }
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_lift; ++i) {
    kkt_residual.lift[i].lx().setRandom();
    kkt_residual.lift[i].lu().setRandom();
    kkt_residual.lift[i].Fx().setRandom();
  }
  return kkt_residual;
}


ContactSequence RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, T, N);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  for (int i=0; i<max_num_impulse; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    tmp.eventTime = i * 0.15 + 0.01 * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
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
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  EXPECT_TRUE(factorization[N].Pqq.isApprox(kkt_matrix[N].Qqq()));
  EXPECT_TRUE(factorization[N].Pqv.isZero());
  EXPECT_TRUE(factorization[N].Pvq.isZero());
  EXPECT_TRUE(factorization[N].Pvv.isApprox(kkt_matrix[N].Qvv()));
  EXPECT_TRUE(factorization[N].sq.isApprox(-1*kkt_residual[N].lq()));
  EXPECT_TRUE(factorization[N].sv.isApprox(-1*kkt_residual[N].lv()));
  const SplitRiccatiFactorization riccati_factorization_default = SplitRiccatiFactorization(robot);
  for (int i=N-1; i>=0; --i) {
    EXPECT_TRUE(factorization[i].isApprox(riccati_factorization_default));
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.impulse[i].isApprox(riccati_factorization_default));
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.aux[i].isApprox(riccati_factorization_default));
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.lift[i].isApprox(riccati_factorization_default));
  }
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  auto factorizer_ref = factorizer;
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  for (int i=N-1; i>=0; --i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref.aux[impulse_index].backwardRiccatiRecursion(
          factorization_ref[i+1], 
          dtau_aux, 
          kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], 
          factorization_ref.aux[impulse_index]);
      factorizer_ref.impulse[impulse_index].backwardRiccatiRecursion(
          factorization_ref.aux[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], 
          factorization_ref.impulse[impulse_index]);
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref.impulse[impulse_index], 
          dtau_impulse, 
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref[i]);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau;
      const double dtau_aux = dtau - dtau_lift;
      ASSERT_TRUE(dtau_lift > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref.lift[lift_index].backwardRiccatiRecursion(
          factorization_ref[i+1], 
          dtau_aux, 
          kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], 
          factorization_ref.lift[lift_index]);
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref.lift[lift_index], 
          dtau_lift, 
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref[i]);
    }
    else {
      factorizer_ref[i].backwardRiccatiRecursion(
          factorization_ref[i+1], 
          dtau, 
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref[i]);
    }
  }
  testIsSame(factorization, factorization_ref);
  testIsSame(kkt_matrix, kkt_matrix_ref);
  testIsSame(kkt_residual, kkt_residual_ref);
}


void RiccatiRecursionTest::testForwardRiccatiRecursionParallel(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  auto constraint_factorization_ref = constraint_factorization;
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, contact_sequence, kkt_matrix, kkt_residual, constraint_factorization);
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  for (int i=0; i<N; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      constraint_factorization_ref.Eq(impulse_index) = kkt_matrix_ref.impulse[impulse_index].Pq();
      constraint_factorization_ref.e(impulse_index) = kkt_residual_ref.impulse[impulse_index].P();
      constraint_factorization_ref.T_impulse(impulse_index, impulse_index).topRows(robot.dimv())
          = kkt_matrix_ref.impulse[impulse_index].Pq().transpose();
      constraint_factorization_ref.T_impulse(impulse_index, impulse_index).bottomRows(robot.dimv()).setZero();
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          exist_state_constraint);
      factorizer_ref.aux[impulse_index].forwardRiccatiRecursionParallel(
          kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], 
          exist_state_constraint);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau;
      const double dtau_aux = dtau - dtau_lift;
      ASSERT_TRUE(dtau_lift > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardRiccatiRecursionParallel(
          kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], 
          exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardRiccatiRecursionParallel(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          exist_state_constraint);
    }
  }
  testIsSame(factorization, factorization_ref);
  testIsSame(kkt_matrix, kkt_matrix_ref);
  testIsSame(kkt_residual, kkt_residual_ref);
  testIsConstraintFactorizationSame(constraint_factorization, constraint_factorization_ref);
}


void RiccatiRecursionTest::testForwardStateConstraintFactorization(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.forwardStateConstraintFactorization(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  for (int i=0; i<N; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i],
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          dtau_impulse,
          factorization_ref.impulse[impulse_index], 
          exist_state_constraint);
      factorizer_ref.impulse[impulse_index].forwardStateConstraintFactorization(
          factorization_ref.impulse[impulse_index], 
          kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], 
          factorization_ref.aux[impulse_index],
          exist_state_constraint);
      factorizer_ref.aux[impulse_index].forwardStateConstraintFactorization(
          factorization_ref.aux[impulse_index],
          kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], 
          dtau_aux,
          factorization_ref[i+1], 
          exist_state_constraint);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau;
      const double dtau_aux = dtau - dtau_lift;
      ASSERT_TRUE(dtau_lift > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i],
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          dtau_lift,
          factorization_ref.lift[lift_index], 
          exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardStateConstraintFactorization(
          factorization_ref.lift[lift_index],
          kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], 
          dtau_aux,
          factorization_ref[i+1], 
          exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardStateConstraintFactorization(
          factorization_ref[i],
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          dtau,
          factorization_ref[i+1], 
          exist_state_constraint);
    }
  }
  testIsSame(factorization, factorization_ref);
  testIsSame(kkt_matrix, kkt_matrix_ref);
  testIsSame(kkt_residual, kkt_residual_ref);
}


void RiccatiRecursionTest::testBackwardStateConstraintFactorization(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, contact_sequence, kkt_matrix, kkt_residual, constraint_factorization);
  riccati_recursion.forwardStateConstraintFactorization(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  auto constraint_factorization_ref = constraint_factorization;
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  riccati_recursion.backwardStateConstraintFactorization(factorizer, contact_sequence, kkt_matrix, constraint_factorization);
  const int num_constraint = contact_sequence.totalNumImpulseStages();
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
    const int time_stage_before_constraint = contact_sequence.timeStageBeforeImpulse(constraint_index);
    const double dtau_constraint = contact_sequence.impulseTime(constraint_index) - time_stage_before_constraint * dtau;
    ASSERT_TRUE(dtau_constraint > 0);
    ASSERT_TRUE(dtau_constraint < dtau);
    factorizer_ref[time_stage_before_constraint].backwardStateConstraintFactorization(
        constraint_factorization_ref.T_impulse(constraint_index, constraint_index), 
        kkt_matrix[time_stage_before_constraint], 
        dtau_constraint, 
        constraint_factorization_ref.T(constraint_index, time_stage_before_constraint));
    for (int i=time_stage_before_constraint-1; i>=0; --i) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau;
        const double dtau_aux = dtau - dtau_impulse;
        ASSERT_TRUE(dtau_impulse > 0);
        ASSERT_TRUE(dtau_aux > 0);
        factorizer_ref.aux[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1),
            kkt_matrix_ref.aux[impulse_index], 
            dtau_aux, 
            constraint_factorization_ref.T_aux(constraint_index, impulse_index));
        factorizer_ref.impulse[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_aux(constraint_index, impulse_index),
            kkt_matrix_ref.impulse[impulse_index], 
            constraint_factorization_ref.T_impulse(constraint_index, impulse_index));
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_impulse(constraint_index, impulse_index),
            kkt_matrix_ref[i], 
            dtau_impulse, 
            constraint_factorization_ref.T(constraint_index, i));
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau;
        const double dtau_aux = dtau - dtau_lift;
        ASSERT_TRUE(dtau_lift > 0);
        ASSERT_TRUE(dtau_aux > 0);
        factorizer_ref.lift[lift_index].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1), 
            kkt_matrix_ref.lift[lift_index], 
            dtau_aux, 
            constraint_factorization_ref.T_lift(constraint_index, lift_index));
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T_lift(constraint_index, lift_index), 
            kkt_matrix_ref[i], 
            dtau_lift, 
            constraint_factorization_ref.T(constraint_index, i));
      }
      else {
        factorizer_ref[i].backwardStateConstraintFactorization(
            constraint_factorization_ref.T(constraint_index, i+1), 
            kkt_matrix_ref[i], 
            dtau, 
            constraint_factorization_ref.T(constraint_index, i));
      }
    }
  }
  testIsSame(factorization, factorization_ref);
  testIsSame(kkt_matrix, kkt_matrix_ref);
  testIsSame(kkt_residual, kkt_residual_ref);
  testIsConstraintFactorizationSame(constraint_factorization, constraint_factorization_ref);
}


void RiccatiRecursionTest::testForwardRiccatiRecursion(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_matrix = createKKTMatrix(robot, contact_sequence);
  auto kkt_residual= createKKTResidual(robot, contact_sequence);
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  RiccatiRecursion riccati_recursion(robot, T, N, max_num_impulse, nproc);
  riccati_recursion.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, factorization);
  RiccatiFactorizer factorizer(N, max_num_impulse, robot);
  riccati_recursion.backwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  riccati_recursion.forwardRiccatiRecursionParallel(factorizer, contact_sequence, kkt_matrix, kkt_residual, constraint_factorization);
  riccati_recursion.forwardStateConstraintFactorization(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization);
  riccati_recursion.backwardStateConstraintFactorization(factorizer, contact_sequence, kkt_matrix, constraint_factorization);
  auto factorizer_ref = factorizer;
  auto factorization_ref = factorization;
  auto kkt_matrix_ref = kkt_matrix; 
  auto kkt_residual_ref = kkt_residual; 
  Direction d(N, max_num_impulse, robot);
  for (int i=0; i<=N; ++i) {
    d[i].dx().setRandom();
  }
  for (int i=0; i<max_num_impulse; ++i) {
    d.impulse[i].dx().setRandom();
    d.aux[i].dx().setRandom();
    d.lift[i].dx().setRandom();
  }
  auto d_ref = d;
  riccati_recursion.forwardRiccatiRecursion(factorizer, contact_sequence, kkt_matrix, kkt_residual, factorization, d);
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  for (int i=0; i<N; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau;
      const double dtau_aux = dtau - dtau_impulse;
      ASSERT_TRUE(dtau_impulse > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref.impulse[impulse_index], 
          d_ref[i],
          dtau_impulse,
          d_ref.impulse[impulse_index],
          exist_state_constraint);
      factorizer_ref.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.impulse[impulse_index], 
          kkt_residual_ref.impulse[impulse_index], 
          d_ref.impulse[impulse_index],
          d_ref.aux[impulse_index]);
      factorizer_ref.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix_ref.aux[impulse_index], 
          kkt_residual_ref.aux[impulse_index], 
          factorization_ref[i+1], 
          d_ref.aux[impulse_index],
          dtau_aux,
          d_ref[i+1],
          exist_state_constraint);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau;
      const double dtau_aux = dtau - dtau_lift;
      ASSERT_TRUE(dtau_lift > 0);
      ASSERT_TRUE(dtau_aux > 0);
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref.lift[lift_index], 
          d_ref[i],
          dtau_lift,
          d_ref.lift[lift_index],
          exist_state_constraint);
      factorizer_ref.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix_ref.lift[lift_index], 
          kkt_residual_ref.lift[lift_index], 
          factorization_ref[i+1], 
          d_ref.lift[lift_index],
          dtau_aux,
          d_ref[i+1],
          exist_state_constraint);
    }
    else {
      factorizer_ref[i].forwardRiccatiRecursion(
          kkt_matrix_ref[i], 
          kkt_residual_ref[i], 
          factorization_ref[i+1], 
          d_ref[i],
          dtau,
          d_ref[i+1],
          exist_state_constraint);
    }
  }
  testIsApprox(d, d_ref);
  testIsSame(factorization, factorization_ref);
  testIsSame(kkt_matrix, kkt_matrix_ref);
  testIsSame(kkt_residual, kkt_residual_ref);
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