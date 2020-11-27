#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

class StateConstraintRiccatiFactorizerTest : public ::testing::Test {
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
  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nproc;
  double T, dtau;
};


ContactSequence StateConstraintRiccatiFactorizerTest::createContactSequence(const Robot& robot) const {
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


void StateConstraintRiccatiFactorizerTest::test(const Robot& robot) const {
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
  std::vector<RiccatiFactorization> impulse_riccati_factorization(max_num_impulse, RiccatiFactorization(robot));
  for (int i=0; i<num_impulse; ++i) {
    const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(dimx, dimx);
    impulse_riccati_factorization[i].N = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    impulse_riccati_factorization[i].Pi.setRandom();
    impulse_riccati_factorization[i].pi.setRandom();
  }
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(dimx);
  std::vector<ImpulseSplitDirection> d(max_num_impulse, ImpulseSplitDirection(robot));
  for (int i=0; i<num_impulse; ++i) {
    d[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d[i].dxi().setRandom();
  }
  StateConstraintRiccatiFactorizer factorizer(robot, max_num_impulse, nproc);
  auto impulse_riccati_factorization_ref = impulse_riccati_factorization;
  auto constraint_factorization_ref = constraint_factorization;
  auto d_ref = d;
  factorizer.computeLagrangeMultiplierDirection(contact_sequence, impulse_riccati_factorization, 
                                                constraint_factorization, dx0, d);
  StateConstraintRiccatiLPFactorizer lp_factorizer(robot);
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    lp_factorizer.factorizeLinearProblem(contact_sequence, impulse_riccati_factorization_ref[constraint_index], 
                                         constraint_factorization_ref, dx0, constraint_index);
  }
  constraint_factorization_ref.ENT().triangularView<Eigen::StrictlyLower>() 
      = constraint_factorization_ref.ENT().transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(constraint_factorization_ref.ENT().isApprox(constraint_factorization.ENT()));
  EXPECT_TRUE(constraint_factorization_ref.e().isApprox(constraint_factorization.e()));
  constraint_factorization_ref.dxi() = constraint_factorization_ref.ENT().inverse() * constraint_factorization_ref.e();
  EXPECT_TRUE(constraint_factorization_ref.dxi().isApprox(constraint_factorization.dxi()));
  for (int i=0; i<num_impulse; ++i) {
    EXPECT_TRUE(d[i].dxi().isApprox(constraint_factorization_ref.dxi(i)));
  }
}


TEST_F(StateConstraintRiccatiFactorizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(StateConstraintRiccatiFactorizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}