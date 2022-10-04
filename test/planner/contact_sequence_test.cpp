#include <vector>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/planner/discrete_event.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/robot/robot.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ContactSequenceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    max_num_each_events = 20;
  }

  virtual void TearDown() {
  }

  static DiscreteEvent createDiscreteEvent(const Robot& robot, 
                                           const ContactStatus& pre_contact_status);
  static std::vector<DiscreteEvent> createDiscreteEvents(const Robot& robot, 
                                                         const ContactStatus& initial_contact_status, 
                                                         const int num_discrete_events);
  void test_constructor(const Robot& robot) const;
  void test_setContactStatus(const Robot& robot) const;
  void test_push_back(const Robot& robot) const;
  void test_pop_back(const Robot& robot) const;
  void test_pop_front(const Robot& robot) const;
  void test_setContactPlacements(const Robot& robot) const;

  int max_num_each_events;
};


DiscreteEvent ContactSequenceTest::createDiscreteEvent(const Robot& robot, 
                                                       const ContactStatus& pre_contact_status) {
  ContactStatus post_contact_status = pre_contact_status;
  post_contact_status.setRandom();
  DiscreteEvent discrete_event(post_contact_status, pre_contact_status);
  while (!discrete_event.existDiscreteEvent()) {
    post_contact_status.setRandom();
    discrete_event.setDiscreteEvent(pre_contact_status, post_contact_status);
  }
  return discrete_event;
}


std::vector<DiscreteEvent> ContactSequenceTest::createDiscreteEvents(const Robot& robot, 
                                                                     const ContactStatus& initial_contact_status, 
                                                                     const int num_discrete_events) {
  std::vector<DiscreteEvent> discrete_event;
  ContactStatus pre_contact_status = initial_contact_status;
  ContactStatus post_contact_status = robot.createContactStatus();
  for (int i=0; i<num_discrete_events; ++i) {
    DiscreteEvent tmp(pre_contact_status, post_contact_status);
    while (!tmp.existDiscreteEvent()) {
      post_contact_status.setRandom();
      tmp.setDiscreteEvent(pre_contact_status, post_contact_status);
    }
    discrete_event.push_back(tmp);
    pre_contact_status = post_contact_status;
  }
  return discrete_event;
}


void ContactSequenceTest::test_constructor(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  auto contact_status = robot.createContactStatus();
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  EXPECT_EQ(contact_sequence.reservedNumDiscreteEvents(), max_num_each_events);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  contact_sequence.pop_back();
  contact_sequence.pop_front();
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  EXPECT_EQ(contact_sequence.reservedNumDiscreteEvents(), max_num_each_events);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
}


void ContactSequenceTest::test_setContactStatus(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  auto default_contact_status = robot.createContactStatus();
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  contact_sequence.init(contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_FALSE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.pop_back();
  EXPECT_FALSE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.init(contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_FALSE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
  contact_sequence.pop_front();
  EXPECT_FALSE(contact_sequence.contactStatus(0) == contact_status);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == default_contact_status);
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_EQ(contact_sequence.numDiscreteEvents(), 0);
  EXPECT_EQ(contact_sequence.numContactPhases(), 1);
}


void ContactSequenceTest::test_push_back(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  auto pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.init(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int j=0; j<5; ++j) {
    contact_sequence.push_back(discrete_events[j], event_times[j], false);
    for (int i=0; i<j+1; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(j+1) == discrete_events[j].postContactStatus());
    int num_impact = 0;
    int num_lift = 0;
    for (int i=0; i<j+1; ++i) {
      EXPECT_EQ(contact_sequence.eventType(num_impact+num_lift), discrete_events[num_impact+num_lift].eventType());
      if (discrete_events[i].existImpact()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impactTime(num_impact), event_times[num_impact+num_lift]);
        EXPECT_TRUE(contact_sequence.impactStatus(num_impact) == discrete_events[num_impact+num_lift].impactStatus());
        ++num_impact;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(num_lift), event_times[num_impact+num_lift]);
        ++num_lift;
      }
    }
    EXPECT_EQ(contact_sequence.numImpactEvents(), num_impact);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
    EXPECT_EQ(contact_sequence.numDiscreteEvents(), num_impact+num_lift);
    EXPECT_EQ(contact_sequence.numContactPhases(), num_impact+num_lift+1);
  }
  EXPECT_NO_THROW(
    std::cout << contact_sequence << std::endl;
  );
  auto contact_sequence_ptr = std::make_shared<ContactSequence>(contact_sequence);
  EXPECT_NO_THROW(
    std::cout << contact_sequence_ptr << std::endl;
  );
}


void ContactSequenceTest::test_pop_back(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.init(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    contact_sequence.push_back(discrete_events[i], event_times[i], false);
  }
  int num_impact = contact_sequence.numImpactEvents();
  int num_lift = contact_sequence.numLiftEvents();
  for (int j=0; j<5; ++j) {
    for (int i=0; i<5-j; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(5-j) == discrete_events[4-j].postContactStatus());
    int num_impact_tmp = 0;
    int num_lift_tmp = 0;
    for (int i=0; i<5-j; ++i) {
      if (discrete_events[i].existImpact()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impactTime(num_impact_tmp), event_times[num_impact_tmp+num_lift_tmp]);
        EXPECT_TRUE(contact_sequence.impactStatus(num_impact_tmp) == discrete_events[num_impact_tmp+num_lift_tmp].impactStatus());
        EXPECT_FALSE(contact_sequence.isSTOEnabledImpact(num_impact_tmp));
        ++num_impact_tmp;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(num_lift_tmp), event_times[num_impact_tmp+num_lift_tmp]);
        EXPECT_FALSE(contact_sequence.isSTOEnabledLift(num_lift_tmp));
        ++num_lift_tmp;
      }
    }
    contact_sequence.pop_back();
    if (discrete_events[5-j-1].existImpact()) {
      --num_impact;
    }
    else {
      --num_lift;
    }
    EXPECT_EQ(contact_sequence.numImpactEvents(), num_impact);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
  }
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
}


void ContactSequenceTest::test_pop_front(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  ContactStatus pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.init(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  for (int i=0; i<5; ++i) {
    contact_sequence.push_back(discrete_events[i], event_times[i], false);
  }
  int num_impact = contact_sequence.numImpactEvents();
  int num_lift = contact_sequence.numLiftEvents();
  for (int j=0; j<5; ++j) {
    for (int i=0; i<5-j; ++i) {
      EXPECT_TRUE(contact_sequence.contactStatus(i) == discrete_events[i+j].preContactStatus());
    }
    EXPECT_TRUE(contact_sequence.contactStatus(5-j) == discrete_events[4].postContactStatus());
    int impact_index = contact_sequence.numImpactEvents() - 1;
    int lift_index = contact_sequence.numLiftEvents() - 1;
    for (int i=4-j; i>=0; --i) {
      if (discrete_events[i+j].existImpact()) {
        EXPECT_DOUBLE_EQ(contact_sequence.impactTime(impact_index), event_times[i+j]);
        EXPECT_TRUE(contact_sequence.impactStatus(impact_index) == discrete_events[i+j].impactStatus());
        --impact_index;
      } 
      else {
        EXPECT_DOUBLE_EQ(contact_sequence.liftTime(lift_index), event_times[i+j]);
        --lift_index;
      }
    }
    contact_sequence.pop_front();
    if (discrete_events[j].existImpact()) {
      --num_impact;
    }
    else {
      --num_lift;
    }
    EXPECT_EQ(contact_sequence.numImpactEvents(), num_impact);
    EXPECT_EQ(contact_sequence.numLiftEvents(), num_lift);
  }
  EXPECT_EQ(contact_sequence.numImpactEvents(), 0);
  EXPECT_EQ(contact_sequence.numLiftEvents(), 0);
  EXPECT_TRUE(contact_sequence.contactStatus(0) == discrete_events[4].postContactStatus());
}


void ContactSequenceTest::test_setContactPlacements(const Robot& robot) const {
  ContactSequence contact_sequence(robot, max_num_each_events);
  auto pre_contact_status = robot.createContactStatus();
  pre_contact_status.setRandom();
  contact_sequence.init(pre_contact_status);
  std::vector<DiscreteEvent> discrete_events = createDiscreteEvents(robot, pre_contact_status, 5);
  std::vector<double> event_times = {0.1, 0.25, 0.5, 0.7, 0.9};
  std::vector<int> impact_indices;
  int impact_index = 0;
  for (int i=0; i<5; ++i) {
    contact_sequence.push_back(discrete_events[i], event_times[i], false);
    if (discrete_events[i].eventType() == DiscreteEventType::Impact) {
      impact_indices.push_back(impact_index);
      ++impact_index;
    }
    else {
      impact_indices.push_back(-1);
    }
  }
  std::vector<std::vector<Eigen::Vector3d>> contact_positions;
  for (int i=0; i<6; ++i) {
    std::vector<Eigen::Vector3d> cp;
    for (int j=0; j<robot.maxNumContacts(); ++j) {
      cp.push_back(Eigen::Vector3d::Random());
    }
    contact_positions.push_back(cp);
  }
  for (int i=0; i<6; ++i) {
    contact_sequence.setContactPlacements(i, contact_positions[i]);
    const auto& cps = contact_sequence.contactStatus(i).contactPositions();
    for (int j=0; j<robot.maxNumContacts(); ++j) {
      EXPECT_TRUE(cps[j].isApprox(contact_positions[i][j]));
    }
    if (i > 0) {
      if (contact_sequence.eventType(i-1) == DiscreteEventType::Impact) {
        const auto& ips = contact_sequence.impactStatus(impact_indices[i-1]).contactPositions();
        for (int j=0; j<robot.maxNumContacts(); ++j) {
          EXPECT_TRUE(ips[j].isApprox(contact_positions[i][j]));
        }
      }
    }
  }
}


TEST_F(ContactSequenceTest, fixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateRobotManipulator(dt);
  test_constructor(robot);
  test_setContactStatus(robot);
  test_push_back(robot);
  test_pop_back(robot);
  test_pop_front(robot);
  test_setContactPlacements(robot);
}


TEST_F(ContactSequenceTest, floatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  test_constructor(robot);
  test_setContactStatus(robot);
  test_push_back(robot);
  test_pop_back(robot);
  test_pop_front(robot);
  test_setContactPlacements(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}