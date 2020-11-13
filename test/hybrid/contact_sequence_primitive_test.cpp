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
  static void testNoDiscreteEvent(const ContactSequencePrimitive& contact_sequence, const int time_stage);
  static void testDiscreteEvent(const ContactSequencePrimitive& contact_sequence, const int time_stage, const DiscreteEvent& discrete_event);
  void testConstructor(const Robot& robot) const;
  void testSetContactStatus(const Robot& robot) const;
  void testSetDiscreteEvent(const Robot& robot) const;
  void testShiftDiscreteEventBeyondInitial(const Robot& robot) const;
  void testShiftDiscreteEventBeyondTerminal(const Robot& robot) const;
  void testShiftDiscreteEventForward(const Robot& robot) const;
  void testShiftDiscreteEventBackward(const Robot& robot) const;
  void testShiftDiscreteEventForwardSameStage(const Robot& robot) const;
  void testShiftDiscreteEventBackwardSameStage(const Robot& robot) const;

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


void ContactSequencePrimitiveTest::testNoDiscreteEvent(const ContactSequencePrimitive& contact_sequence, const int time_stage) {
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(time_stage), 0);
  EXPECT_FALSE(contact_sequence.existDiscreteEvent(time_stage));
  EXPECT_FALSE(contact_sequence.existImpulse(time_stage));
  EXPECT_FALSE(contact_sequence.existLift(time_stage));
  EXPECT_FALSE(contact_sequence.existOnlyLift(time_stage));
}


void ContactSequencePrimitiveTest::testDiscreteEvent(const ContactSequencePrimitive& contact_sequence, const int time_stage, const DiscreteEvent& discrete_event) {
  EXPECT_TRUE(contact_sequence.contactStatus(time_stage) == discrete_event.preContactStatus());
  EXPECT_TRUE(contact_sequence.impulseStatus(time_stage) == discrete_event.impulseStatus());
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(time_stage), 0);
  EXPECT_EQ(contact_sequence.existDiscreteEvent(time_stage), discrete_event.existDiscreteEvent());
  EXPECT_EQ(contact_sequence.existImpulse(time_stage), discrete_event.existImpulse());
  EXPECT_EQ(contact_sequence.existLift(time_stage), discrete_event.existLift());
  EXPECT_EQ(contact_sequence.existOnlyLift(time_stage), (discrete_event.existLift() && (!discrete_event.existImpulse())));
}


void ContactSequencePrimitiveTest::testConstructor(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus contact_status = robot.createContactStatus();
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=1; i<5; ++i) {
    for (int j=event_time_stages[i-1]+1; j<event_time_stages[i]; ++j) {
      EXPECT_TRUE(contact_sequence.contactStatus(j) == discrete_events[i].preContactStatus());
      testNoDiscreteEvent(contact_sequence, j);
    }
    EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[i]));
    testDiscreteEvent(contact_sequence, event_time_stages[i], discrete_events[i]);
  }
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
  discrete_events[1] = createDiscreteEvent(robot, discrete_events[0].postContactStatus());
  const int new_event_time_stage = 6;
  contact_sequence.setDiscreteEvent(discrete_events[1], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(new_event_time_stage));
  testDiscreteEvent(contact_sequence, new_event_time_stage, discrete_events[1]);
  for (int i=new_event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[3]));
  testDiscreteEvent(contact_sequence, event_time_stages[3], discrete_events[3]);
  for (int i=event_time_stages[3]+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  testDiscreteEvent(contact_sequence, event_time_stages[4], discrete_events[4]);
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<event_time_stages[1]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[1]));
  testDiscreteEvent(contact_sequence, event_time_stages[1], discrete_events[1]);
  for (int i=event_time_stages[1]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(new_event_time_stage) == discrete_events[0].postContactStatus());
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[0].postContactStatus(), discrete_events[3].postContactStatus());
  testDiscreteEvent(contact_sequence, new_event_time_stage, new_discrete_event);
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  testDiscreteEvent(contact_sequence, event_time_stages[4], discrete_events[4]);
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
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
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[1].preContactStatus(), discrete_events[4].preContactStatus());
  testDiscreteEvent(contact_sequence, new_event_time_stage, new_discrete_event);
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  testDiscreteEvent(contact_sequence, event_time_stages[4], discrete_events[4]);
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventForwardSameStage(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  const int new_event_time_stage = event_time_stages[1];
  contact_sequence.shiftDiscreteEvent(event_time_stages[3], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(new_event_time_stage) == discrete_events[0].postContactStatus());
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[0].postContactStatus(), discrete_events[4].preContactStatus());
  testDiscreteEvent(contact_sequence, new_event_time_stage, new_discrete_event);
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  testDiscreteEvent(contact_sequence, event_time_stages[4], discrete_events[4]);
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.contactStatus(N) == discrete_events[4].postContactStatus());
}


void ContactSequencePrimitiveTest::testShiftDiscreteEventBackwardSameStage(const Robot& robot) const {
  ContactSequencePrimitive contact_sequence(robot, N);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const std::vector<int> event_time_stages = {4, 7, 8, 9, 15};
  for (int i=0; i<5; ++i) {
    contact_sequence.setDiscreteEvent(discrete_events[i], event_time_stages[i]);
  }
  const int new_event_time_stage = event_time_stages[3];
  contact_sequence.shiftDiscreteEvent(event_time_stages[1], new_event_time_stage);
  for (int i=0; i<event_time_stages[0]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[0].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[0]));
  testDiscreteEvent(contact_sequence, event_time_stages[0], discrete_events[0]);
  for (int i=event_time_stages[0]+1; i<new_event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[1].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  DiscreteEvent new_discrete_event(robot);
  new_discrete_event.setDiscreteEvent(discrete_events[0].postContactStatus(), discrete_events[4].preContactStatus());
  testDiscreteEvent(contact_sequence, new_event_time_stage, new_discrete_event);
  for (int i=new_event_time_stage+1; i<event_time_stages[4]; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].preContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
  }
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stages[4]));
  testDiscreteEvent(contact_sequence, event_time_stages[4], discrete_events[4]);
  for (int i=event_time_stages[4]+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[4].postContactStatus());
    testNoDiscreteEvent(contact_sequence, i);
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
  testShiftDiscreteEventForwardSameStage(robot);
  testShiftDiscreteEventBackwardSameStage(robot);
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
  testShiftDiscreteEventForwardSameStage(robot);
  testShiftDiscreteEventBackwardSameStage(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}