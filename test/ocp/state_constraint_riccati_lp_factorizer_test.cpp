#include <string>
#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

class StateConstraintRiccatiLPFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    T = 1;
    N = 20;
    max_num_impulse = 5;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;
  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  double T, dtau;
  int N, max_num_impulse;
};


ContactSequence StateConstraintRiccatiLPFactorizerTest::createContactSequence(const Robot& robot) const {
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
    tmp.eventTime = i * 0.15 + 0.1 * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_impulse; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void StateConstraintRiccatiLPFactorizerTest::test(const Robot& robot) const {
  const int dimv = robot.dimv();
  const auto contact_sequence = createContactSequence(robot);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  std::vector<RiccatiFactorization> impulse_riccati_factorization(max_num_impulse, RiccatiFactorization(robot));
  for (int i=0; i<num_impulse; ++i) {
    const int dimx = 2*robot.dimv();
    const Eigen::MatrixXd seed_mat = Eigen::MatrixXd::Random(dimx, dimx);
    impulse_riccati_factorization[i].N = seed_mat * seed_mat.transpose() + Eigen::MatrixXd::Identity(dimx, dimx);
    impulse_riccati_factorization[i].Pi.setRandom();
    impulse_riccati_factorization[i].pi.setRandom();
  }
  StateConstraintRiccatiFactorization constraint_factorization(robot, N, max_num_impulse);
  constraint_factorization.setConstraintStatus(contact_sequence);
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    for (int i=0; i<max_num_impulse; ++i) {
      constraint_factorization.T_impulse(constraint_index, i).setRandom();
    }
    constraint_factorization.Eq(constraint_index).setRandom();
    constraint_factorization.e(constraint_index).setRandom();
    constraint_factorization.T_impulse(constraint_index, constraint_index).topRows(dimv)
        = constraint_factorization.Eq(constraint_index).transpose();
    constraint_factorization.T_impulse(constraint_index, constraint_index).bottomRows(dimv).setZero();
  }
  auto constraint_factorization_ref = constraint_factorization;
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*robot.dimv());
  StateConstraintRiccatiLPFactorizer factorizer(robot);
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    factorizer.factorizeLinearProblem(contact_sequence, impulse_riccati_factorization[constraint_index], 
                                      constraint_factorization, dx0, constraint_index);
  }
  for (int constraint_index=0; constraint_index<num_impulse; ++constraint_index) {
    const auto impulse_riccati_factorization_ref = impulse_riccati_factorization[constraint_index];
    const int dimf = contact_sequence.impulseStatus(constraint_index).dimp();
    const int dimv = robot.dimv();
    const int dimx = 2*robot.dimv();
    Eigen::MatrixXd E = Eigen::MatrixXd::Zero(dimf, dimx);
    E.leftCols(dimv) = constraint_factorization_ref.Eq(constraint_index);
    ASSERT_TRUE(impulse_riccati_factorization_ref.N.isApprox(impulse_riccati_factorization_ref.N.transpose()));
    const Eigen::MatrixXd EN = E * impulse_riccati_factorization_ref.N;
    const Eigen::MatrixXd ENEt = EN * E.transpose();
    const Eigen::VectorXd e = constraint_factorization_ref.e(constraint_index) 
                              + E * impulse_riccati_factorization_ref.Pi * dx0
                              + E * impulse_riccati_factorization_ref.pi;
    EXPECT_TRUE(EN.isApprox(constraint_factorization.EN(constraint_index)));
    EXPECT_TRUE(ENEt.isApprox(constraint_factorization.ENEt(constraint_index)));
    EXPECT_TRUE(e.isApprox(constraint_factorization.e(constraint_index)));
    EXPECT_TRUE(ENEt.llt().info() == Eigen::Success);
    EXPECT_TRUE(constraint_factorization.ENEt(constraint_index).isApprox(constraint_factorization.ENEt(constraint_index).transpose()));
    constraint_factorization_ref.EN(constraint_index) = EN;
    constraint_factorization_ref.ENEt(constraint_index) = ENEt;
    constraint_factorization_ref.e(constraint_index) = e;
    for (int i=constraint_index+1; i<num_impulse; ++i) {
      constraint_factorization_ref.ENT(constraint_index, i) 
          = EN * constraint_factorization_ref.T_impulse(i, constraint_index);
    }
  }
  EXPECT_TRUE(constraint_factorization_ref.ENT().isApprox(constraint_factorization.ENT()));
  EXPECT_TRUE(constraint_factorization_ref.e().isApprox(constraint_factorization.e()));
}


TEST_F(StateConstraintRiccatiLPFactorizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(StateConstraintRiccatiLPFactorizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}