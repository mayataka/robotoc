#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/hybrid/discrete_event.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/robot/robot.hpp"


namespace idocp {

class ContactSequenceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    T = 1;
    dtau = T / N;
  }

  virtual void TearDown() {
  }

  static DiscreteEvent createDiscreteEvent(const Robot& robot, const ContactStatus& pre_contact_status);
  static std::vector<DiscreteEvent> createDiscreteEvents(const Robot& robot, const ContactStatus& initial_contact_status, const int num);
  void testConstructor(const Robot& robot) const;
  void testSetContactStatus(const Robot& robot) const;
  void testSetDiscreteEvent(const Robot& robot) const;
  void testShiftDiscreteEventBeyondInitial(const Robot& robot) const;
  void testShiftDiscreteEventBeyondTerminal(const Robot& robot) const;
  void testShiftDiscreteEventInitial(const Robot& robot) const;
  void testShiftDiscreteEventTerminal(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N;
  double T, dtau;
};


DiscreteEvent ContactSequenceTest::createDiscreteEvent(const Robot& robot, const ContactStatus& pre_contact_status) {
  DiscreteEvent discrete_event(robot);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  while (!discrete_event.existDiscreteEvent()) {
    post_contact_status.setRandom();
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


std::vector<DiscreteEvent> ContactSequenceTest::createDiscreteEvents(const Robot& robot, const ContactStatus& initial_contact_status, const int num) {
  std::vector<DiscreteEvent> discrete_event;
  ContactStatus pre_contact_status = initial_contact_status;
  ContactStatus post_contact_status = robot.createContactStatus();
  std::random_device rnd;
  for (int i=0; i<num; ++i) {
    DiscreteEvent tmp(robot);
    tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    discrete_event.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  return discrete_event;
}


void ContactSequenceTest::testConstructor(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == contact_status);
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


void ContactSequenceTest::testSetContactStatus(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  std::random_device rnd;
  contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(contact_status);
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == contact_status);
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


void ContactSequenceTest::testSetDiscreteEvent(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    discrete_events[i].eventTime = event_times[i];
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  int num_impulse = 0;
  int num_lift = 0;
  for (int i=0; i<5; ++i) {
    const int stage = std::floor(event_times[i] / dtau);
    EXPECT_EQ(contact_sequence.numImpulseStages(stage), num_impulse);
    EXPECT_EQ(contact_sequence.numLiftStages(stage), num_lift);
    if (discrete_events[i].existImpulse()) {
      EXPECT_EQ(contact_sequence.numImpulseStages(stage+1), num_impulse+1);
      EXPECT_EQ(contact_sequence.numLiftStages(stage+1), num_lift);
      EXPECT_EQ(contact_sequence.timeStageBeforeImpulse(num_impulse), stage);
      EXPECT_DOUBLE_EQ(contact_sequence.impulseTime(num_impulse), event_times[num_impulse+num_lift]);
      EXPECT_TRUE(contact_sequence.impulseStatus(num_impulse) == discrete_events[num_impulse+num_lift].impulseStatus());
      ++num_impulse;
    } 
    else {
      EXPECT_EQ(contact_sequence.numImpulseStages(stage+1), num_impulse);
      EXPECT_EQ(contact_sequence.numLiftStages(stage+1), num_lift+1);
      EXPECT_EQ(contact_sequence.timeStageBeforeLift(num_lift), stage);
      EXPECT_DOUBLE_EQ(contact_sequence.liftTime(num_lift), event_times[num_impulse+num_lift]);
      ++num_lift;
    }
  }
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), num_impulse);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), num_lift);
}


void ContactSequenceTest::testShiftDiscreteEventBeyondInitial(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    discrete_events[i].eventTime = event_times[i];
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  if (contact_sequence.totalNumImpulseStages() > 0) {
    contact_sequence.shiftImpulse(contact_sequence.totalNumImpulseStages()-1, -1);
    EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  }
  if (contact_sequence.totalNumLiftStages() > 0) {
    contact_sequence.shiftLift(contact_sequence.totalNumLiftStages()-1, -1);
    EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  }
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


void ContactSequenceTest::testShiftDiscreteEventBeyondTerminal(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    discrete_events[i].eventTime = event_times[i];
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  if (contact_sequence.totalNumImpulseStages() > 0) {
    contact_sequence.shiftImpulse(0, T+1);
    EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  }
  if (contact_sequence.totalNumLiftStages() > 0) {
    contact_sequence.shiftLift(0, T+1);
    EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  }
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[0].preContactStatus());
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


void ContactSequenceTest::testShiftDiscreteEventInitial(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    discrete_events[i].eventTime = event_times[i];
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  if (contact_sequence.totalNumImpulseStages() > 0) {
    contact_sequence.shiftImpulse(contact_sequence.totalNumImpulseStages()-1, 0);
    EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  }
  if (contact_sequence.totalNumLiftStages() > 0) {
    contact_sequence.shiftLift(contact_sequence.totalNumLiftStages()-1, 0);
    EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  }
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


void ContactSequenceTest::testShiftDiscreteEventTerminal(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    discrete_events[i].eventTime = event_times[i];
    contact_sequence.setDiscreteEvent(discrete_events[i]);
  }
  if (contact_sequence.totalNumImpulseStages() > 0) {
    contact_sequence.shiftImpulse(0, T);
    EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  }
  if (contact_sequence.totalNumLiftStages() > 0) {
    contact_sequence.shiftLift(0, T);
    EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  }
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.numImpulseStages(i), 0);
    EXPECT_EQ(contact_sequence.numLiftStages(i), 0);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[0].preContactStatus());
  EXPECT_EQ(contact_sequence.totalNumImpulseStages(), 0);
  EXPECT_EQ(contact_sequence.totalNumLiftStages(), 0);
}


TEST_F(ContactSequenceTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testSetDiscreteEvent(robot);
  testShiftDiscreteEventBeyondInitial(robot);
  testShiftDiscreteEventBeyondTerminal(robot);
  testShiftDiscreteEventInitial(robot);
  testShiftDiscreteEventTerminal(robot);
}


TEST_F(ContactSequenceTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testSetDiscreteEvent(robot);
  testShiftDiscreteEventBeyondInitial(robot);
  testShiftDiscreteEventBeyondTerminal(robot);
  testShiftDiscreteEventInitial(robot);
  testShiftDiscreteEventTerminal(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}