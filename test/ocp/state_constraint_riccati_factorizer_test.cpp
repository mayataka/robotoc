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
    N = 10;
    max_num_impulse = N;
    T = 1;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  static ContactStatus createRamdomContactStatus(const Robot& robot);
  void testComputeLagrangeMultiplierDirection(const Robot& robot) const;
  void testAggregateLagrangeMultiplierDirection(const Robot& robot) const;
  void testAggregateLagrangeMultiplierDirectionImpulse(const Robot& robot) const;
  void testAggregateLagrangeMultiplierDirectionLift(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse;
  double T, dtau;
};


ContactStatus StateConstraintRiccatiFactorizerTest::createRamdomContactStatus(const Robot& robot) {
  ContactStatus contact_status = robot.createContactStatus();
  contact_status.setRandom();
  return contact_status;
}


void StateConstraintRiccatiFactorizerTest::testComputeLagrangeMultiplierDirection(const Robot& robot) const {
  std::vector<RiccatiFactorization> impulse_riccati_factorization(max_num_impulse, RiccatiFactorization(robot));
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization(max_num_impulse, StateConstraintRiccatiFactorization(robot, N, max_num_impulse));
  std::vector<ImpulseSplitDirection> d(max_num_impulse, ImpulseSplitDirection(robot));
  std::vector<StateConstraintRiccatiLPFactorizer> lp_factorizer(max_num_impulse, StateConstraintRiccatiLPFactorizer(robot));
  ContactSequence contact_sequence(robot, T, N);
  const int num_proc = 4;
  StateConstraintRiccatiFactorizer factorizer(robot, max_num_impulse, num_proc);
  ContactStatus pre_contact_status = robot.createContactStatus();
  for (int i=0; i<max_num_impulse; ++i) {
    ContactStatus post_contact_status = createRamdomContactStatus(robot);
    DiscreteEvent discrete_event(robot);
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
    discrete_event.eventTime = T*i/max_num_impulse + 0.5*dtau;
    if (discrete_event.existDiscreteEvent()) {
      contact_sequence.setDiscreteEvent(discrete_event);
    }
    pre_contact_status = post_contact_status;
  }
  std::cout << "num impulse = " << contact_sequence.totalNumImpulseStages() << std::endl;
  for (int i=0; i<contact_sequence.totalNumImpulseStages(); ++i) {
    impulse_riccati_factorization[i].N.setRandom();
    impulse_riccati_factorization[i].N = (impulse_riccati_factorization[i].N * impulse_riccati_factorization[i].N.transpose()).eval();
    impulse_riccati_factorization[i].N += Eigen::MatrixXd::Identity(2*robot.dimv(), 2*robot.dimv());
    impulse_riccati_factorization[i].Pi.setRandom();
    impulse_riccati_factorization[i].pi.setRandom();
    constraint_factorization[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    constraint_factorization[i].Eq().setRandom();
    constraint_factorization[i].e().setRandom();
    d[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    for (int j=0; j<contact_sequence.totalNumImpulseStages(); ++j) {
      constraint_factorization[i].T_impulse(j).setRandom();
    }
  }
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization_ref = constraint_factorization;
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*robot.dimv());
  factorizer.computeLagrangeMultiplierDirection(contact_sequence, impulse_riccati_factorization, 
                                                constraint_factorization, dx0, d);
  int dim_total = 0;
  std::vector<int> dims;
  for (int i=0; i<contact_sequence.totalNumImpulseStages(); ++i) {
    lp_factorizer[i].factorizeLinearProblem(impulse_riccati_factorization[i], constraint_factorization_ref[i], dx0);
    const int dim = contact_sequence.impulseStatus(i).dimp();
    dims.push_back(dim);
    dim_total += dim;
  }
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(dim_total, dim_total);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(dim_total);
  int block_begin = 0;
  for (int i=0; i<contact_sequence.totalNumImpulseStages(); ++i) {
    const int dim = dims[i];
    A.block(block_begin, block_begin, dim, dim) = constraint_factorization_ref[i].ENEt();
    int col_begin = block_begin + dim;
    for (int j=i+1; j<contact_sequence.totalNumImpulseStages(); ++j) {
      const int dim_col = dims[j];
      A.block(block_begin, col_begin, dim, dim_col) 
          = constraint_factorization_ref[i].EN() * constraint_factorization_ref[j].T_impulse(i);
      col_begin += dim_col;
    }
    b.segment(block_begin, dim) = constraint_factorization_ref[i].e();
    block_begin += dim;
  }
  const Eigen::VectorXd dxi_ref = A.inverse() * b;
  block_begin = 0;
  for (int i=0; i<contact_sequence.totalNumImpulseStages(); ++i) {
    const int dim = dims[i];
    EXPECT_TRUE(d[i].dxi().isApprox(dxi_ref.segment(block_begin, dim)));
    block_begin += dim;
  }
}


void StateConstraintRiccatiFactorizerTest::testAggregateLagrangeMultiplierDirection(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  const int num_proc = 4;
  ContactStatus pre_contact_status = robot.createContactStatus();
  for (int i=0; i<max_num_impulse; ++i) {
    ContactStatus post_contact_status = createRamdomContactStatus(robot);
    DiscreteEvent discrete_event(robot);
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
    discrete_event.eventTime = T*i/max_num_impulse + 0.5*dtau;
    if (discrete_event.existDiscreteEvent()) {
      contact_sequence.setDiscreteEvent(discrete_event);
    }
    pre_contact_status = post_contact_status;
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  std::vector<ImpulseSplitDirection> d(max_num_impulse, ImpulseSplitDirection(robot));
  for (int i=0; i<num_impulse; ++i) {
    d[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d[i].dxi().setRandom();
  }
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization(max_num_impulse, StateConstraintRiccatiFactorization(robot, N, max_num_impulse));
  for (int i=0; i<num_impulse; ++i) {
    constraint_factorization[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    for (int j=0; j<N; ++j) {
      constraint_factorization[i].T(j).setRandom();
    }
  }
  std::vector<RiccatiFactorization> riccati_factorization(N, RiccatiFactorization(robot));
  for (int i=0; i<N; ++i) {
    StateConstraintRiccatiFactorizer::aggregateLagrangeMultiplierDirection(
        contact_sequence, constraint_factorization, d, i, riccati_factorization[i]);
  }
  std::vector<RiccatiFactorization> riccati_factorization_ref(N, RiccatiFactorization(robot));
  for (int i=0; i<N; ++i) {
    riccati_factorization_ref[i].n.setZero();
    for (int j=num_impulse-1; j>=0; --j) {
      const int time_stage_before_impulse = contact_sequence.timeStageBeforeImpulse(j);
      if (time_stage_before_impulse >= i) {
        riccati_factorization_ref[i].n += constraint_factorization[j].T(i) * d[j].dxi();
      }
      else {
        break;
      }
    }
  }
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(riccati_factorization[i].n.isApprox(riccati_factorization_ref[i].n));
  }
}


void StateConstraintRiccatiFactorizerTest::testAggregateLagrangeMultiplierDirectionImpulse(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  const int num_proc = 4;
  ContactStatus pre_contact_status = robot.createContactStatus();
  for (int i=0; i<max_num_impulse; ++i) {
    ContactStatus post_contact_status = createRamdomContactStatus(robot);
    DiscreteEvent discrete_event(robot);
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
    discrete_event.eventTime = T*i/max_num_impulse + 0.5*dtau;
    if (discrete_event.existDiscreteEvent()) {
      contact_sequence.setDiscreteEvent(discrete_event);
    }
    pre_contact_status = post_contact_status;
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  std::vector<ImpulseSplitDirection> d(max_num_impulse, ImpulseSplitDirection(robot));
  for (int i=0; i<num_impulse; ++i) {
    d[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d[i].dxi().setRandom();
  }
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization(max_num_impulse, StateConstraintRiccatiFactorization(robot, N, max_num_impulse));
  for (int i=0; i<num_impulse; ++i) {
    constraint_factorization[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    for (int j=0; j<num_impulse; ++j) {
      constraint_factorization[i].T_impulse(j).setRandom();
    }
  }
  std::vector<RiccatiFactorization> riccati_factorization(max_num_impulse, RiccatiFactorization(robot));
  for (int i=0; i<num_impulse; ++i) {
    StateConstraintRiccatiFactorizer::aggregateLagrangeMultiplierDirectionImpulse(
        contact_sequence, constraint_factorization, d, i, riccati_factorization[i]);
  }
  std::vector<RiccatiFactorization> riccati_factorization_ref(max_num_impulse, RiccatiFactorization(robot));
  for (int i=0; i<max_num_impulse; ++i) {
    riccati_factorization_ref[i].n.setZero();
    for (int j=num_impulse-1; j>=i; --j) {
      riccati_factorization_ref[i].n += constraint_factorization[j].T_impulse(i) * d[j].dxi();
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    EXPECT_TRUE(riccati_factorization[i].n.isApprox(riccati_factorization_ref[i].n));
  }
}


void StateConstraintRiccatiFactorizerTest::testAggregateLagrangeMultiplierDirectionLift(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  const int num_proc = 4;
  ContactStatus pre_contact_status = robot.createContactStatus();
  for (int i=0; i<max_num_impulse; ++i) {
    ContactStatus post_contact_status = createRamdomContactStatus(robot);
    DiscreteEvent discrete_event(robot);
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
    discrete_event.eventTime = T*i/max_num_impulse + 0.5*dtau;
    if (discrete_event.existDiscreteEvent()) {
      contact_sequence.setDiscreteEvent(discrete_event);
    }
    pre_contact_status = post_contact_status;
  }
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  std::vector<ImpulseSplitDirection> d(max_num_impulse, ImpulseSplitDirection(robot));
  for (int i=0; i<num_impulse; ++i) {
    d[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    d[i].dxi().setRandom();
  }
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization(max_num_impulse, StateConstraintRiccatiFactorization(robot, N, max_num_impulse));
  const int max_num_lift = max_num_impulse;
  const int num_lift = contact_sequence.totalNumLiftStages();
  for (int i=0; i<num_impulse; ++i) {
    constraint_factorization[i].setImpulseStatus(contact_sequence.impulseStatus(i));
    for (int j=0; j<num_lift; ++j) {
      constraint_factorization[i].T_lift(j).setRandom();
    }
  }
  std::vector<RiccatiFactorization> riccati_factorization(max_num_lift, RiccatiFactorization(robot));
  for (int i=0; i<num_lift; ++i) {
    StateConstraintRiccatiFactorizer::aggregateLagrangeMultiplierDirectionLift(
        contact_sequence, constraint_factorization, d, i, riccati_factorization[i]);
  }
  std::vector<RiccatiFactorization> riccati_factorization_ref(max_num_lift, RiccatiFactorization(robot));
  for (int i=0; i<num_lift; ++i) {
    const int time_stage_before_lift = contact_sequence.timeStageBeforeLift(i);
    riccati_factorization_ref[i].n.setZero();
    for (int j=num_impulse-1; j>=0; --j) {
      const int time_stage_before_impulse = contact_sequence.timeStageBeforeImpulse(j);
      if (time_stage_before_impulse >= time_stage_before_lift) {
        riccati_factorization_ref[i].n += constraint_factorization[j].T_lift(i) * d[j].dxi();
      }
      else {
        break;
      }
    }
  }
  for (int i=0; i<max_num_lift; ++i) {
    EXPECT_TRUE(riccati_factorization[i].n.isApprox(riccati_factorization_ref[i].n));
  }
}


TEST_F(StateConstraintRiccatiFactorizerTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testComputeLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirectionImpulse(robot);
  testAggregateLagrangeMultiplierDirectionLift(robot);
}


TEST_F(StateConstraintRiccatiFactorizerTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testComputeLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirection(robot);
  testAggregateLagrangeMultiplierDirectionImpulse(robot);
  testAggregateLagrangeMultiplierDirectionLift(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}