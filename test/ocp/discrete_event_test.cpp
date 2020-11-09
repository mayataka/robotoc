#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/discrete_event.hpp"


namespace idocp {

class DiscreteEventTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    max_point_contacts = 10;
  }

  virtual void TearDown() {
  }

  int max_point_contacts;
};


TEST_F(DiscreteEventTest, constructor) {
  DiscreteEvent discrete_event(max_point_contacts);
  EXPECT_EQ(discrete_event.max_point_contacts(), max_point_contacts);
  EXPECT_FALSE(discrete_event.hasDiscreteEvent());
  EXPECT_FALSE(discrete_event.hasImpulse());
  EXPECT_FALSE(discrete_event.hasLift());
  EXPECT_DOUBLE_EQ(discrete_event.eventTime(), 0);
  ContactStatus contact_status(max_point_contacts);
  EXPECT_EQ(contact_status.num_active_contacts(), 0);
  discrete_event.act(contact_status);
  EXPECT_EQ(contact_status.num_active_contacts(), 0);
  const double event_time = 10;
  discrete_event.setEventTime(event_time);
  EXPECT_DOUBLE_EQ(event_time, discrete_event.eventTime());
}


TEST_F(DiscreteEventTest, impulse) {
  DiscreteEvent discrete_event(max_point_contacts);
  ContactStatus cs_before(max_point_contacts), cs_after(max_point_contacts);
  cs_before.activateContacts({1, 2, 3});
  cs_after.activateContacts({1, 2, 3, 4, 5, 6});
  discrete_event.setDiscreteEvent(cs_before, cs_after);
  EXPECT_EQ(discrete_event.max_point_contacts(), max_point_contacts);
  EXPECT_TRUE(discrete_event.hasDiscreteEvent());
  EXPECT_TRUE(discrete_event.hasImpulse());
  EXPECT_FALSE(discrete_event.hasLift());
  EXPECT_DOUBLE_EQ(discrete_event.eventTime(), 0);
  discrete_event.act(cs_before);
  EXPECT_TRUE(cs_before == cs_after);
  discrete_event.disableDiscreteEvent();
  EXPECT_FALSE(discrete_event.hasDiscreteEvent());
  EXPECT_FALSE(discrete_event.hasImpulse());
  EXPECT_FALSE(discrete_event.hasLift());
}


TEST_F(DiscreteEventTest, lift) {
  DiscreteEvent discrete_event(max_point_contacts);
  ContactStatus cs_before(max_point_contacts), cs_after(max_point_contacts);
  cs_before.activateContacts({1, 2, 3, 4, 5, 6});
  cs_after.activateContacts({1, 2, 3});
  discrete_event.setDiscreteEvent(cs_before, cs_after);
  EXPECT_EQ(discrete_event.max_point_contacts(), max_point_contacts);
  EXPECT_TRUE(discrete_event.hasDiscreteEvent());
  EXPECT_FALSE(discrete_event.hasImpulse());
  EXPECT_TRUE(discrete_event.hasLift());
  EXPECT_DOUBLE_EQ(discrete_event.eventTime(), 0);
  discrete_event.act(cs_before);
  EXPECT_TRUE(cs_before == cs_after);
  discrete_event.disableDiscreteEvent();
  EXPECT_FALSE(discrete_event.hasDiscreteEvent());
  EXPECT_FALSE(discrete_event.hasImpulse());
  EXPECT_FALSE(discrete_event.hasLift());
}


TEST_F(DiscreteEventTest, impulseAndLift) {
  DiscreteEvent discrete_event(max_point_contacts);
  ContactStatus cs_before(max_point_contacts), cs_after(max_point_contacts);
  cs_before.activateContacts({1, 2, 3, 6, 7});
  cs_after.activateContacts({3, 5, 6, 7, 8, 9});
  discrete_event.setDiscreteEvent(cs_before, cs_after);
  EXPECT_EQ(discrete_event.max_point_contacts(), max_point_contacts);
  EXPECT_TRUE(discrete_event.hasDiscreteEvent());
  EXPECT_TRUE(discrete_event.hasImpulse());
  EXPECT_TRUE(discrete_event.hasLift());
  EXPECT_DOUBLE_EQ(discrete_event.eventTime(), 0);
  discrete_event.act(cs_before);
  EXPECT_TRUE(cs_before == cs_after);
  discrete_event.disableDiscreteEvent();
  EXPECT_FALSE(discrete_event.hasDiscreteEvent());
  EXPECT_FALSE(discrete_event.hasImpulse());
  EXPECT_FALSE(discrete_event.hasLift());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}