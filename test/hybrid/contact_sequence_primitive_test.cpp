#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/hybrid/discrete_event.hpp"
#include "idocp/hybrid/contact_sequence_primitive.hpp"
#include "idocp/robot/robot.hpp"


namespace idocp {

class ContactSequencePrimitiveTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
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
  void testShiftDiscreteEventForward(const Robot& robot) const;
  void testShiftDiscreteEventBackward(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N;
};


DiscreteEvent ContactSequencePrimitiveTest::createDiscreteEvent(const Robot& robot, const ContactStatus& pre_contact_status) {
  DiscreteEvent discrete_event(robot);
  ContactStatus post_contact_status = pre_contact_status;
  std::random_device rnd;
  while (!discrete_event.existDiscreteEvent()) {
    post_contact_status.setRandom();
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


std::vector<DiscreteEvent> ContactSequencePrimitiveTest::createDiscreteEvents(const Robot& robot, const ContactStatus& initial_contact_status, const int num) {
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


void ContactSequencePrimitiveTest::testConstructor(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == contact_status);
}


void ContactSequencePrimitiveTest::testSetContactStatus(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  std::random_device rnd;
  contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(contact_status);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == contact_status);
}


void ContactSequencePrimitiveTest::testSetDiscreteEvent(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[0]) == discrete_events[0].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[0]) == discrete_events[0].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[0]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[0]), discrete_events[0].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[0]), discrete_events[0].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[0]), (discrete_events[0].existLift() && (!discrete_events[0].existImpulse())));
  for (int i=1; i<5; ++i) {
    for (int j=event_time_stages[i-1]+1; j<event_time_stages[i]; ++j) {
      EXPECT_TRUE(contact_sequence.contactStatus(j) == discrete_events[i].preContactStatus());
      EXPECT_EQ(contact_sequence.eventTime(j), 0);
      EXPECT_FALSE(contact_sequence.existDiscreteEvent(j));
      EXPECT_FALSE(contact_sequence.existImpulse(j));
      EXPECT_FALSE(contact_sequence.existLift(j));
      EXPECT_FALSE(contact_sequence.existOnlyLift(j));
    }
    EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[i]) == discrete_events[i].preContactStatus());
    EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[i]) == discrete_events[i].impulseStatus());
    EXPECT_EQ(contact_sequence.eventTime(event_time_stages[i]), 0);
    EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[i]));
    EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[i]), discrete_events[i].existImpulse());
    EXPECT_EQ(contact_sequence.existLift(event_time_stages[i]), discrete_events[i].existLift());
    EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[i]), (discrete_events[i].existLift() && (!discrete_events[i].existImpulse())));
  }
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
  discrete_events[1] = createDiscreteEvent(robot, discrete_events[0].postContactStatus());
  const int new_event_time_stage = 6;
  contact_sequence.setDiscreteEvent(discrete_events[1], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[0]) == discrete_events[0].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[0]) == discrete_events[0].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[0]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[0]), discrete_events[0].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[0]), discrete_events[0].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[0]), (discrete_events[0].existLift() && (!discrete_events[0].existImpulse())));
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(new_event_time_stage) == discrete_events[1].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(new_event_time_stage) == discrete_events[1].impulseStatus());
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(new_event_time_stage), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(new_event_time_stage));
  EXPECT_EQ(contact_sequence.existImpulse(new_event_time_stage), discrete_events[1].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(new_event_time_stage), discrete_events[1].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(new_event_time_stage), (discrete_events[1].existLift() && (!discrete_events[1].existImpulse())));
  for (int i=new_event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[1].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventBeyondInitial(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  contact_sequence.shiftDiscreteEventBeyondInitial(event_time_stages[2]);
  for (int i=0; i<event_time_stages[3]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[3].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[3]) == discrete_events[3].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[3]) == discrete_events[3].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[3]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[3]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[3]), discrete_events[3].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[3]), discrete_events[3].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[3]), (discrete_events[3].existLift() && (!discrete_events[3].existImpulse())));
  for (int i=event_time_stages[3]+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[4]) == discrete_events[4].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[4]) == discrete_events[4].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[4]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[4]), discrete_events[4].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[4]), discrete_events[4].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[4]), (discrete_events[4].existLift() && (!discrete_events[4].existImpulse())));
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventBeyondTerminal(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  contact_sequence.shiftDiscreteEventBeyondTerminal(event_time_stages[2]);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[0]) == discrete_events[0].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[0]) == discrete_events[0].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[0]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[0]), discrete_events[0].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[0]), discrete_events[0].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[0]), (discrete_events[0].existLift() && (!discrete_events[0].existImpulse())));
  for (int i=event_time_stages[0]+1; i<event_time_stages[1]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[1]) == discrete_events[1].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[1]) == discrete_events[1].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[1]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[1]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[1]), discrete_events[1].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[1]), discrete_events[1].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[1]), (discrete_events[1].existLift() && (!discrete_events[1].existImpulse())));
  for (int i=event_time_stages[1]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[1].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventForward(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  const int new_event_time_stage = 6;
  contact_sequence.shiftDiscreteEvent(event_time_stages[3], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[0]) == discrete_events[0].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[0]) == discrete_events[0].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[0]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[0]), discrete_events[0].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[0]), discrete_events[0].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[0]), (discrete_events[0].existLift() && (!discrete_events[0].existImpulse())));
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(new_event_time_stage) == discrete_events[0].postContactStatus());
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[0].postContactStatus(), discrete_events[3].postContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(new_event_time_stage) == new_discrete_event.impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(new_event_time_stage), 0);
  EXPECT_EQ(contact_sequence.existDiscreteEvent(new_event_time_stage), new_discrete_event.existDiscreteEvent());
  EXPECT_EQ(contact_sequence.existImpulse(new_event_time_stage), new_discrete_event.existImpulse());
  EXPECT_EQ(contact_sequence.existLift(new_event_time_stage), new_discrete_event.existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(new_event_time_stage), (new_discrete_event.existLift() && (!new_discrete_event.existImpulse())));
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[4]) == discrete_events[4].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[4]) == discrete_events[4].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[4]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[4]), discrete_events[4].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[4]), discrete_events[4].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[4]), (discrete_events[4].existLift() && (!discrete_events[4].existImpulse())));
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventBackward(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  const int new_event_time_stage = 14;
  contact_sequence.shiftDiscreteEvent(event_time_stages[1], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[0]) == discrete_events[0].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[0]) == discrete_events[0].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[0]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[0]), discrete_events[0].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[0]), discrete_events[0].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[0]), (discrete_events[0].existLift() && (!discrete_events[0].existImpulse())));
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[1].preContactStatus(), discrete_events[4].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(new_event_time_stage) == new_discrete_event.impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(new_event_time_stage), 0);
  EXPECT_EQ(contact_sequence.existDiscreteEvent(new_event_time_stage), new_discrete_event.existDiscreteEvent());
  EXPECT_EQ(contact_sequence.existImpulse(new_event_time_stage), new_discrete_event.existImpulse());
  EXPECT_EQ(contact_sequence.existLift(new_event_time_stage), new_discrete_event.existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(new_event_time_stage), (new_discrete_event.existLift() && (!new_discrete_event.existImpulse())));
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stages[4]) == discrete_events[4].preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(event_time_stages[4]) == discrete_events[4].impulseStatus());
  EXPECT_EQ(contact_sequence.eventTime(event_time_stages[4]), 0);
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  EXPECT_EQ(contact_sequence.existImpulse(event_time_stages[4]), discrete_events[4].existImpulse());
  EXPECT_EQ(contact_sequence.existLift(event_time_stages[4]), discrete_events[4].existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(event_time_stages[4]), (discrete_events[4].existLift() && (!discrete_events[4].existImpulse())));
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    EXPECT_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
    EXPECT_FALSE(contact_sequence.existOnlyLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
}


TEST_F(ContactSequencePrimitiveTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testSetDiscreteEvent(robot);
  testShiftDiscreteEventBeyondInitial(robot);
  testShiftDiscreteEventBeyondTerminal(robot);
  testShiftDiscreteEventForward(robot);
  testShiftDiscreteEventBackward(robot);
}


TEST_F(ContactSequencePrimitiveTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testSetDiscreteEvent(robot);
  testShiftDiscreteEventBeyondInitial(robot);
  testShiftDiscreteEventBeyondTerminal(robot);
  testShiftDiscreteEventForward(robot);
  testShiftDiscreteEventBackward(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}