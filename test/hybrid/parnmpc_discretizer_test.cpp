#include <random>
#include <vector>
#include <algorithm>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"
#include "idocp/robot/robot.hpp"

namespace idocp {

class ParNMPCDiscretizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_events = 5;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    T = 1;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  ContactSequence createContactSequence(const Robot& robot) const;

  void testConstructor(const Robot& robot) const;
  void testDiscretizeOCP(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_events;
  double t, T, dtau;
};


ContactSequence ParNMPCDiscretizerTest::createContactSequence(const Robot& robot) const {
  std::vector<DiscreteEvent> discrete_events;
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  ContactSequence contact_sequence(robot, N);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  for (int i=0; i<max_num_events; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    const double event_period = 3 * dtau;
    tmp.eventTime = t + i * event_period + 0.1 * event_period * std::abs(Eigen::VectorXd::Random(1)[0]);
    discrete_events.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  for (int i=0; i<max_num_events; ++i) {
    contact_sequence.pushBackDiscreteEvent(discrete_events[i]);
  }
  return contact_sequence;
}


void ParNMPCDiscretizerTest::testConstructor(const Robot& robot) const {
  ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_events);
  EXPECT_EQ(parnmpc_discretizer.N(), N);
  EXPECT_EQ(parnmpc_discretizer.numImpulseStages(), 0);
  EXPECT_EQ(parnmpc_discretizer.numLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(parnmpc_discretizer.contactPhase(i), 0);
  }
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageBeforeImpulse(i));
    EXPECT_FALSE(parnmpc_discretizer.isTimeStageBeforeLift(i));
  }
}

void ParNMPCDiscretizerTest::testDiscretizeOCP(const Robot& robot) const {
  ParNMPCDiscretizer parnmpc_discretizer(T, N, max_num_events);
  OCPDiscretizer ocp_discretizer(T, N, max_num_events);
  const ContactSequence contact_sequence = createContactSequence(robot);
  parnmpc_discretizer.discretizeOCP(contact_sequence, t);
  ocp_discretizer.discretizeOCP(contact_sequence, t);
  EXPECT_EQ(parnmpc_discretizer.N(), ocp_discretizer.N());
  EXPECT_EQ(parnmpc_discretizer.numImpulseStages(), ocp_discretizer.numImpulseStages());
  EXPECT_EQ(parnmpc_discretizer.numLiftStages(), ocp_discretizer.numLiftStages());
  for (int i=0; i<N; ++i) {
    EXPECT_EQ(parnmpc_discretizer.contactPhase(i), ocp_discretizer.contactPhase(i+1));
  }
  for (int i=0; i<N; ++i) {
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.t(i), ocp_discretizer.t(i+1));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.dtau(i), ocp_discretizer.dtau(i+1));
    if (parnmpc_discretizer.isTimeStageBeforeImpulse(i)) {
      EXPECT_TRUE(parnmpc_discretizer.isTimeStageAfterImpulse(i+1));
      EXPECT_TRUE(ocp_discretizer.isTimeStageBeforeImpulse(i+1));
      EXPECT_EQ(parnmpc_discretizer.impulseIndexAfterTimeStage(i), ocp_discretizer.impulseIndexAfterTimeStage(i+1));
    }
    if (parnmpc_discretizer.isTimeStageBeforeLift(i)) {
      EXPECT_TRUE(parnmpc_discretizer.isTimeStageAfterLift(i+1));
      EXPECT_TRUE(ocp_discretizer.isTimeStageBeforeLift(i+1));
      EXPECT_EQ(parnmpc_discretizer.liftIndexAfterTimeStage(i), ocp_discretizer.liftIndexAfterTimeStage(i+1));
    }
  }
  for (int i=0; i<contact_sequence.numImpulseEvents(); ++i) {
    EXPECT_EQ(parnmpc_discretizer.timeStageBeforeImpulse(i)+1, ocp_discretizer.timeStageBeforeImpulse(i));
    EXPECT_EQ(parnmpc_discretizer.timeStageAfterImpulse(i)+1, ocp_discretizer.timeStageAfterImpulse(i));
    EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeImpulse(i), ocp_discretizer.contactPhaseBeforeImpulse(i));
    EXPECT_EQ(parnmpc_discretizer.contactPhaseAfterImpulse(i), ocp_discretizer.contactPhaseAfterImpulse(i));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.t_impulse(i), ocp_discretizer.t_impulse(i));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.dtau_aux(i), ocp_discretizer.dtau_aux(i));
  }
  for (int i=0; i<contact_sequence.numLiftEvents(); ++i) {
    EXPECT_EQ(parnmpc_discretizer.timeStageBeforeLift(i)+1, ocp_discretizer.timeStageBeforeLift(i));
    EXPECT_EQ(parnmpc_discretizer.timeStageAfterLift(i)+1, ocp_discretizer.timeStageAfterLift(i));
    EXPECT_EQ(parnmpc_discretizer.contactPhaseBeforeLift(i), ocp_discretizer.contactPhaseBeforeLift(i));
    EXPECT_EQ(parnmpc_discretizer.contactPhaseAfterLift(i), ocp_discretizer.contactPhaseAfterLift(i));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.t_lift(i), ocp_discretizer.t_lift(i));
    EXPECT_DOUBLE_EQ(parnmpc_discretizer.dtau_lift(i), ocp_discretizer.dtau_lift(i));
  }
}


TEST_F(ParNMPCDiscretizerTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testDiscretizeOCP(robot);
}


TEST_F(ParNMPCDiscretizerTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testDiscretizeOCP(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}