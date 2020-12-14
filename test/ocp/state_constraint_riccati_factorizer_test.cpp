#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"


namespace idocp {

class StateConstraintSplitRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_impulse = 5;
    T = 1;
    dtau = T / N;
    nproc = 4;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;
  void testComputeLagrangeMultiplierDirection(const Robot& robot) const;
  void testAggregateLagrangeMultiplierDirection(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, dtau;
};


ContactSequence StateConstraintSplitRiccatiFactorizerTest::createContactSequence(const Robot& robot) const {
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


void StateConstraintSplitRiccatiFactorizerTest::testComputeLagrangeMultiplierDirection(const Robot& robot) const {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const auto contact_sequence = createContactSequence(robot);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  for (int i=0; i<num_impulse; ++i) {
    for (int j=0; j<i; ++j) {
      constraint_factorization.T_impulse(i, j).setRandom();
    }
    constraint_factorization.Eq(i).setRandom();
    constraint_factorization.e(i).setRandom();
    constraint_factorization.T_impulse(i, i).topRows(dimv) = constraint_factorization.Eq(i).transpose();
    constraint_factorization.T_impulse(i, i).bottomRows(dimv).setZero();
  }
  RiccatiFactorization riccati_factorization(N, max_num_impulse, robot);
  for (int i=0; i<num_impulse; ++i) {
    const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(dimx, dimx);
    riccati_factorization.impulse[i].N = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    riccati_factorization.impulse[i].Pi.setRandom();
    riccati_factorization.impulse[i].pi.setRandom();
  }
  Direction d(N, max_num_impulse, robot);
  d[0].dx().setRandom();
  for (int i=0; i<num_impulse; ++i) {
    d.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d.impulse[i].dxi().setRandom();
  }
  StateConstraintRiccatiFactorizer factorizer(robot, N, max_num_impulse, nproc);
  auto riccati_factorization_ref = riccati_factorization;
  auto constraint_factorization_ref = constraint_factorization;
  auto d_ref = d;
  factorizer.computeLagrangeMultiplierDirection(contact_sequence, riccati_factorization, 
                                                constraint_factorization, d);
  StateConstraintRiccatiLPFactorizer lp_factorizer(robot);
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    lp_factorizer.factorizeLinearProblem(contact_sequence, riccati_factorization_ref.impulse[constraint_index], 
                                         constraint_factorization_ref, d_ref[0].dx(), constraint_index);
  }
  constraint_factorization_ref.ENT().triangularView<Eigen::StrictlyLower>() 
      = constraint_factorization_ref.ENT().transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(constraint_factorization_ref.ENT().isApprox(constraint_factorization.ENT()));
  EXPECT_TRUE(constraint_factorization_ref.e().isApprox(constraint_factorization.e()));
  constraint_factorization_ref.dxi() = constraint_factorization_ref.ENT().inverse() * constraint_factorization_ref.e();
  EXPECT_TRUE(constraint_factorization_ref.dxi().isApprox(constraint_factorization.dxi()));
  for (int i=0; i<num_impulse; ++i) {
    EXPECT_TRUE(d.impulse[i].dxi().isApprox(constraint_factorization_ref.dxi(i)));
  }
}


void StateConstraintSplitRiccatiFactorizerTest::testAggregateLagrangeMultiplierDirection(const Robot& robot) const {
  const auto contact_sequence = createContactSequence(robot);
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_impulse; ++i) {
    for (int j=0; j<N; ++j) {
      constraint_factorization.T(i, j).setRandom();
    }
    for (int j=0; j<num_impulse; ++j) {
      constraint_factorization.T_impulse(i, j).setRandom();
      constraint_factorization.T_aux(i, j).setRandom();
    }
    for (int j=0; j<num_lift; ++j) {
      constraint_factorization.T_lift(i, j).setRandom();
    }
  }
  Direction d(N, max_num_impulse, robot);
  for (int i=0; i<num_impulse; ++i) {
    d.impulse[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d.impulse[i].dxi().setRandom();
  }
  RiccatiFactorization factorization(N, max_num_impulse, robot);
  for (int i=0; i<N; ++i) {
    factorization[i].n.setRandom();
  }
  for (int i=0; i<num_impulse; ++i) {
    factorization.impulse[i].n.setRandom();
    factorization.aux[i].n.setRandom();
  }
  for (int i=0; i<num_lift; ++i) {
    factorization.lift[i].n.setRandom();
  }
  StateConstraintRiccatiFactorizer factorizer(robot, N, max_num_impulse, nproc);
  auto factorization_ref = factorization;
  factorizer.aggregateLagrangeMultiplierDirection(
      constraint_factorization, contact_sequence, d, factorization);
  for (int i=0; i<N; ++i) {
    factorization_ref[i].n.setZero();
    for (int j=0; j<num_impulse; ++j) {
      if (i <= contact_sequence.timeStageBeforeImpulse(j)) {
        factorization_ref[i].n += constraint_factorization.T(j, i) * d.impulse[j].dxi();
      }
    }
  }
  for (int i=0; i<num_impulse; ++i) {
    factorization_ref.impulse[i].n.setZero();
    for (int j=i; j<num_impulse; ++j) {
      factorization_ref.impulse[i].n += constraint_factorization.T_impulse(j, i) * d.impulse[j].dxi();
    }
    factorization_ref.aux[i].n.setZero();
    for (int j=i+1; j<num_impulse; ++j) {
      factorization_ref.aux[i].n += constraint_factorization.T_aux(j, i) * d.impulse[j].dxi();
    }
  }
  for (int i=0; i<num_lift; ++i) {
    const int time_stage_after_lift = contact_sequence.timeStageAfterLift(i);
    factorization_ref.lift[i].n.setZero();
    for (int j=0; j<num_impulse; ++j) {
      if (time_stage_after_lift <= contact_sequence.timeStageBeforeImpulse(j)) {
        factorization_ref.lift[i].n += constraint_factorization.T_lift(j, i) * d.impulse[j].dxi();
      }
    }
  }
  for (int i=0; i<=N; ++i) {
    EXPECT_TRUE(factorization[i].isApprox(factorization_ref[i]));
    EXPECT_FALSE(factorization[i].hasNaN());
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.impulse[i].isApprox(factorization_ref.impulse[i]));
    EXPECT_FALSE(factorization.impulse[i].hasNaN());
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.aux[i].isApprox(factorization_ref.aux[i]));
    EXPECT_FALSE(factorization.aux[i].hasNaN());
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(factorization.lift[i].isApprox(factorization_ref.lift[i]));
    EXPECT_FALSE(factorization.lift[i].hasNaN());
  }
}


TEST_F(StateConstraintSplitRiccatiFactorizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testComputeLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirection(robot);
}


TEST_F(StateConstraintSplitRiccatiFactorizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testComputeLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirection(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}