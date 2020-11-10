#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/discrete_event.hpp"
#include "idocp/ocp/contact_sequence.hpp"
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

  static DiscreteEvent createDiscreteEvent(const Robot& robot, ContactStatus& pre_contact_status);
  void testConstructor(const Robot& robot) const;
  void testSetContactStatus(const Robot& robot) const;
  void testSetDiscreteEvent(const Robot& robot) const;
  void testShiftDiscreteEventToInitial1(const Robot& robot) const;
  void testShiftDiscreteEventToInitial2(const Robot& robot) const;
  void testShiftDiscreteEventToTerminal1(const Robot& robot) const;
  void testShiftDiscreteEventToTerminal2(const Robot& robot) const;
  void testShiftDiscreteEventForward(const Robot& robot) const;
  void testShiftDiscreteEventBackward(const Robot& robot) const;
  void testShiftDiscreteEventMultiple(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N;
  double T, dtau;
};


DiscreteEvent ContactSequenceTest::createDiscreteEvent(const Robot& robot, ContactStatus& pre_contact_status) {
  DiscreteEvent discrete_event(robot.max_point_contacts());
  ContactStatus post_contact_status(robot.max_point_contacts());
  std::random_device rnd;
  while (!discrete_event.existDiscreteEvent()) {
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      post_contact_status.deactivateContact(i);
    }
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (rnd()%2==0) {
        pre_contact_status.activateContact(i);
      }
    }
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (rnd()%2==0) {
        post_contact_status.activateContact(i);
      }
    }
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


void ContactSequenceTest::testConstructor(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status(robot.max_point_contacts());
  ImpulseStatus impulse_status(robot.max_point_contacts());
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testSetContactStatus(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status(robot.max_point_contacts());
  ImpulseStatus impulse_status(robot.max_point_contacts());
  std::random_device rnd;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (rnd()%2==0) {
      contact_status.activateContact(i);
    }
  }
  contact_sequence.setContactStatusUniformly(contact_status);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testSetDiscreteEvent(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  for (int i=0; i<event_time_stage; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(event_time_stage)));
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stage));
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(event_time_stage), (T/2.0));
  if (discrete_event.existImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(event_time_stage));
  }
  if (discrete_event.existLift()) {
    EXPECT_TRUE(contact_sequence.existLift(event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(event_time_stage));
  }
  for (int i=event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventToInitial1(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  contact_sequence.shiftDiscreteEventToInitial(event_time_stage);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_FALSE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventToInitial2(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  contact_sequence.shiftDiscreteEvent(event_time_stage, -1.0);
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_FALSE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventToTerminal1(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  contact_sequence.shiftDiscreteEventToTerminal(event_time_stage);
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventToTerminal2(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  contact_sequence.shiftDiscreteEvent(event_time_stage, T+1.0);
  for (int i=0; i<N; ++i) {
    EXPECT_FALSE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventForward(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  const double shifted_time = T / 3.0;
  const int shifted_event_time_stage = std::floor(shifted_time/dtau);
  contact_sequence.shiftDiscreteEvent(event_time_stage, shifted_time);
  for (int i=0; i<shifted_event_time_stage; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(shifted_event_time_stage)));
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(shifted_event_time_stage));
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(shifted_event_time_stage), shifted_time);
  if (discrete_event.existImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(shifted_event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(shifted_event_time_stage));
  }
  if (discrete_event.existLift()) {
    EXPECT_TRUE(contact_sequence.existLift(shifted_event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(shifted_event_time_stage));
  }
  for (int i=shifted_event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventBackward(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event = createDiscreteEvent(robot, pre_contact_status);
  ASSERT_TRUE(discrete_event.isConsisitentWithPreContactStatus(pre_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  discrete_event.setEventTime(T / 2.0);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2.0)/dtau);
  const double shifted_time = 2.0 * T / 3.0;
  const int shifted_event_time_stage = std::floor(shifted_time/dtau);
  contact_sequence.shiftDiscreteEvent(event_time_stage, shifted_time);
  for (int i=0; i<shifted_event_time_stage; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(discrete_event.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(shifted_event_time_stage)));
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(shifted_event_time_stage));
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(shifted_event_time_stage), shifted_time);
  if (discrete_event.existImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(shifted_event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(shifted_event_time_stage));
  }
  if (discrete_event.existLift()) {
    EXPECT_TRUE(contact_sequence.existLift(shifted_event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(shifted_event_time_stage));
  }
  for (int i=shifted_event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(discrete_event.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEventMultiple(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus pre_contact_status(robot.max_point_contacts()), 
                intermediate_contact_status(robot.max_point_contacts()),
                post_contact_status(robot.max_point_contacts());
  DiscreteEvent discrete_event1(robot.max_point_contacts()),
                discrete_event2(robot.max_point_contacts());
  std::random_device rnd;
  while (!discrete_event1.existDiscreteEvent() || !discrete_event2.existDiscreteEvent()) {
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      pre_contact_status.deactivateContact(i);
      intermediate_contact_status.deactivateContact(i);
      post_contact_status.deactivateContact(i);
    }
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (rnd()%2==0) {
        pre_contact_status.activateContact(i);
      }
    }
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (rnd()%2==0) {
        intermediate_contact_status.activateContact(i);
      }
    }
    for (int i=0; i<robot.max_point_contacts(); ++i) {
      if (rnd()%2==0) {
        post_contact_status.activateContact(i);
      }
    }
    discrete_event1.setDiscreteEvent(pre_contact_status, intermediate_contact_status);
    discrete_event2.setDiscreteEvent(intermediate_contact_status, post_contact_status);
  }
  ASSERT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(pre_contact_status));
  ASSERT_TRUE(discrete_event1.isConsisitentWithPostContactStatus(intermediate_contact_status));
  ASSERT_TRUE(discrete_event2.isConsisitentWithPreContactStatus(intermediate_contact_status));
  ASSERT_TRUE(discrete_event2.isConsisitentWithPostContactStatus(post_contact_status));
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  const double time1 = T / 3.0;
  const double time2 = 2.0 * T / 3.0;
  discrete_event1.setEventTime(time1);
  discrete_event2.setEventTime(time2);
  contact_sequence.setDiscreteEvent(discrete_event1);
  contact_sequence.setDiscreteEvent(discrete_event2);
  const int event_time_stage1 = std::floor(time1/dtau);
  const int event_time_stage2 = std::floor(time2/dtau);
  for (int i=0; i<event_time_stage1; ++i) {
    EXPECT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(event_time_stage1)));
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stage1));
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(event_time_stage1), time1);
  if (discrete_event1.existImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(event_time_stage1));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(event_time_stage1));
  }
  if (discrete_event1.existLift()) {
    EXPECT_TRUE(contact_sequence.existLift(event_time_stage1));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(event_time_stage1));
  }
  for (int i=event_time_stage1+1; i<event_time_stage2; ++i) {
    EXPECT_TRUE(discrete_event1.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_TRUE(discrete_event2.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(discrete_event1.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(event_time_stage2)));
  EXPECT_TRUE(discrete_event2.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(event_time_stage2)));
  EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stage2));
  EXPECT_DOUBLE_EQ(contact_sequence.eventTime(event_time_stage2), time2);
  if (discrete_event2.existImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(event_time_stage2));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(event_time_stage2));
  }
  if (discrete_event2.existLift()) {
    EXPECT_TRUE(contact_sequence.existLift(event_time_stage2));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(event_time_stage2));
  }
  for (int i=event_time_stage2+1; i<N; ++i) {
    EXPECT_TRUE(discrete_event2.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  const double shifted_time1 = 3.0 * T / 4.0;
  ASSERT_TRUE(shifted_time1 > time2);
  const int shifted_event_time_stage = std::floor(shifted_time1/dtau);
  contact_sequence.shiftDiscreteEvent(event_time_stage1, shifted_time1);
  for (int i=0; i<shifted_event_time_stage; ++i) {
    EXPECT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  if (pre_contact_status == post_contact_status) {
    EXPECT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(shifted_event_time_stage)));
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(shifted_event_time_stage));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(shifted_event_time_stage), 0);
  }
  else {
    EXPECT_TRUE(discrete_event1.isConsisitentWithPreContactStatus(contact_sequence.contactStatus(shifted_event_time_stage)));
    EXPECT_TRUE(contact_sequence.existDiscreteEvent(shifted_event_time_stage));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(shifted_event_time_stage), shifted_time1);
    DiscreteEvent discrete_event(robot.max_point_contacts());
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
    if (discrete_event.existImpulse()) {
      EXPECT_TRUE(contact_sequence.existImpulse(shifted_event_time_stage));
    }
    else {
      EXPECT_FALSE(contact_sequence.existImpulse(shifted_event_time_stage));
    }
    if (discrete_event.existLift()) {
      EXPECT_TRUE(contact_sequence.existLift(shifted_event_time_stage));
    }
    else {
      EXPECT_FALSE(contact_sequence.existLift(shifted_event_time_stage));
    }
  }
  for (int i=shifted_event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(discrete_event2.isConsisitentWithPostContactStatus(contact_sequence.contactStatus(i)));
    EXPECT_DOUBLE_EQ(contact_sequence.eventTime(i), 0);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}



TEST_F(ContactSequenceTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testShiftDiscreteEventToInitial1(robot);
  testShiftDiscreteEventToInitial2(robot);
  testShiftDiscreteEventToTerminal1(robot);
  testShiftDiscreteEventToTerminal2(robot);
  testShiftDiscreteEventForward(robot);
  testShiftDiscreteEventBackward(robot);
  testShiftDiscreteEventMultiple(robot);
}


TEST_F(ContactSequenceTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
  testShiftDiscreteEventToInitial1(robot);
  testShiftDiscreteEventToInitial2(robot);
  testShiftDiscreteEventToTerminal1(robot);
  testShiftDiscreteEventToTerminal2(robot);
  testShiftDiscreteEventForward(robot);
  testShiftDiscreteEventBackward(robot);
  testShiftDiscreteEventMultiple(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}