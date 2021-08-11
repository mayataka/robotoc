#include <vector>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/hybrid/discrete_event.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/robot/robot.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ContactSequenceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    max_num_events = 20;
  }

  virtual void TearDown() {
  }

  static DiscreteEvent createDiscreteEvent(const Robot& robot, 
                                           const ContactStatus& pre_contact_status);
  static std::vector<DiscreteEvent> createDiscreteEvents(const Robot& robot, 
                                                         const ContactStatus& initial_contact_status, 
                                                         const int num);
  void testConstructor(const Robot& robot) const;
  void testSetContactStatus(const Robot& robot) const;
  void test_push_back(const Robot& robot) const;
  void test_pop_back(const Robot& robot) const;
  void test_pop_front(const Robot& robot) const;

  int max_num_events;
};


DiscreteEvent ContactSequenceTest::createDiscreteEvent(const Robot& robot, 
                                                       const ContactStatus& pre_contact_status) {
  DiscreteEvent discrete_event(robot.maxPointContacts());
  ContactStatus post_contact_status = pre_contact_status;
  while (!discrete_event.existDiscreteEvent()) {
    post_contact_status.setRandom();
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


std::vector<DiscreteEvent> ContactSequenceTest::createDiscreteEvents(const Robot& robot, 
                                                                     const ContactStatus& initial_contact_status, 
                                                                     const int num) {
  std::vector<DiscreteEvent> discrete_event;
  ContactStatus pre_contact_status = initial_contact_status;
  ContactStatus post_contact_status = robot.createContactStatus();
  for (int i=0; i<num; ++i) {
    DiscreteEvent tmp(robot.maxPointContacts());
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
  ContactSequence contact_sequence(robot, max_num_events);
  auto contact_status = robot.createContactStatus();
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  contact_sequence.pop_back();
  contact_sequence.pop_front();
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
}


void ContactSequenceTest::testSetContactStatus(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_events);
  auto default_contact_status = robot.createContactStatus();
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  contact_sequence.setContactStatusUniformly(contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_FALSE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.pop_back();
  EXPECT_FALSE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.setContactStatusUniformly(contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_FALSE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.pop_front();
  EXPECT_FALSE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
}


void ContactSequenceTest::test_push_back(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_events);
  auto pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int j=0; j<5; ++j) {
    contact_sequence.push_back(discrete_events[j], event_times[j]);
    for (int i=0; i<j+1; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(j+1) == discrete_events[j].postContactStatus());
    int num_impulse = 0;
    int num_lift = 0;
    for (int i=0; i<j+1; ++i) {
      EXPECT_EQ(contact_sequence.eventType(num_impulse+num_lift), discrete_events[num_impulse+num_lift].eventType());
      if (discrete_events[i].existImpulse()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impulseTime(num_impulse), event_times[num_impulse+num_lift]);
        EXPECT_TRUE(contact_sequence.impulseStatus(num_impulse) == discrete_events[num_impulse+num_lift].impulseStatus());
        ++num_impulse;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(num_lift), event_times[num_impulse+num_lift]);
        ++num_lift;
      }
    }
    EXPECT_EQ(contact_sequence.numImpulseEvents(), num_impulse);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
    EXPECT_EQ(contact_sequence.numDiscreteEvents(), num_impulse+num_lift);
    EXPECT_EQ(contact_sequence.numContactPhases(), num_impulse+num_lift+1);
  }
}


void ContactSequenceTest::test_pop_back(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_events);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    contact_sequence.push_back(discrete_events[i], event_times[i]);
  }
  int num_impulse = contact_sequence.numImpulseEvents();
  int num_lift = contact_sequence.numLiftEvents();
  for (int j=0; j<5; ++j) {
    for (int i=0; i<5-j; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(5-j) == discrete_events[4-j].postContactStatus());
    int num_impulse_tmp = 0;
    int num_lift_tmp = 0;
    for (int i=0; i<5-j; ++i) {
      if (discrete_events[i].existImpulse()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impulseTime(num_impulse_tmp), event_times[num_impulse_tmp+num_lift_tmp]);
        EXPECT_TRUE(contact_sequence.impulseStatus(num_impulse_tmp) == discrete_events[num_impulse_tmp+num_lift_tmp].impulseStatus());
        ++num_impulse_tmp;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(num_lift_tmp), event_times[num_impulse_tmp+num_lift_tmp]);
        ++num_lift_tmp;
      }
    }
    contact_sequence.pop_back();
    if (discrete_events[5-j-1].existImpulse()) {
      --num_impulse;
    }
    else {
      --num_lift;
    }
    EXPECT_EQ(contact_sequence.numImpulseEvents(), num_impulse);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
  }
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
}


void ContactSequenceTest::test_pop_front(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_events);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.setContactStatusUniformly(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    contact_sequence.push_back(discrete_events[i], event_times[i]);
  }
  int num_impulse = contact_sequence.numImpulseEvents();
  int num_lift = contact_sequence.numLiftEvents();
  for (int j=0; j<5; ++j) {
    for (int i=0; i<5-j; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i+j].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(5-j) == discrete_events[4].postContactStatus());
    int impulse_index = contact_sequence.numImpulseEvents() - 1;
    int lift_index = contact_sequence.numLiftEvents() - 1;
    for (int i=4-j; i>=0; --i) {
      if (discrete_events[i+j].existImpulse()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impulseTime(impulse_index), event_times[i+j]);
        EXPECT_TRUE(contact_sequence.impulseStatus(impulse_index) == discrete_events[i+j].impulseStatus());
        --impulse_index;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(lift_index), event_times[i+j]);
        --lift_index;
      }
    }
    contact_sequence.pop_front();
    if (discrete_events[j].existImpulse()) {
      --num_impulse;
    }
    else {
      --num_lift;
    }
    EXPECT_EQ(contact_sequence.numImpulseEvents(), num_impulse);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
  }
  EXPECT_EQ(contact_sequence.numImpulseEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == discrete_events[4].postContactStatus());
}


TEST_F(ContactSequenceTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testConstructor(robot);
  testSetContactStatus(robot);
  test_push_back(robot);
  test_pop_back(robot);
  test_pop_front(robot);
}


TEST_F(ContactSequenceTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testConstructor(robot);
  testSetContactStatus(robot);
  test_push_back(robot);
  test_pop_back(robot);
  test_pop_front(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}