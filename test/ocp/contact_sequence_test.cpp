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

  void testConstructor(const Robot& robot) const;
  void testSetContactStatus(const Robot& robot) const;
  void testSetDiscreteEvent(const Robot& robot) const;
  void testShiftDiscreteEvent1(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N;
  double T, dtau;
};


void ContactSequenceTest::testConstructor(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status(robot.max_point_contacts());
  ImpulseStatus impulse_status(robot.max_point_contacts());
  for (int i=0; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status);
    EXPECT_TRUE(contact_sequence.impulseStatus(i) == impulse_status);
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
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testSetDiscreteEvent(const Robot& robot) const {
  ContactSequence contact_sequence(robot, T, N);
  ContactStatus contact_status_before(robot.max_point_contacts());
  ContactStatus contact_status_after(robot.max_point_contacts());
  std::random_device rnd;
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (rnd()%2==0) {
      contact_status_before.activateContact(i);
    }
  }
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    if (rnd()%2==0) {
      contact_status_after.activateContact(i);
    }
  }
  contact_sequence.setContactStatusUniformly(contact_status_before);
  DiscreteEvent discrete_event(robot.max_point_contacts());
  discrete_event.setDiscreteEvent(contact_status_before, contact_status_after);
  discrete_event.setEventTime(T / 2);
  contact_sequence.setDiscreteEvent(discrete_event);
  const int event_time_stage = std::floor((T/2)/dtau);
  for (int i=0; i<event_time_stage; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status_before);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
  EXPECT_TRUE(contact_sequence.contactStatus(event_time_stage) == contact_status_before);
  if (discrete_event.hasDiscreteEvent()) {
    EXPECT_TRUE(contact_sequence.existDiscreteEvent(event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(event_time_stage));
  }
  if (discrete_event.hasImpulse()) {
    EXPECT_TRUE(contact_sequence.existImpulse(event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existImpulse(event_time_stage));
  }
  if (discrete_event.hasLift()) {
    EXPECT_TRUE(contact_sequence.existLift(event_time_stage));
  }
  else {
    EXPECT_FALSE(contact_sequence.existLift(event_time_stage));
  }
  for (int i=event_time_stage+1; i<N; ++i) {
    EXPECT_TRUE(contact_sequence.contactStatus(i) == contact_status_after);
    EXPECT_FALSE(contact_sequence.existDiscreteEvent(i));
    EXPECT_FALSE(contact_sequence.existImpulse(i));
    EXPECT_FALSE(contact_sequence.existLift(i));
  }
}


void ContactSequenceTest::testShiftDiscreteEvent1(const Robot& robot) const {
}


TEST_F(ContactSequenceTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
}


TEST_F(ContactSequenceTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  testConstructor(robot);
  testSetContactStatus(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}