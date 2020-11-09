#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/robot/discrete_event.hpp"


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
  DiscreteEvent impulse_status(max_point_contacts);
  ContactStatus contact_status(max_point_contacts);
  EXPECT_EQ(contact_status.max_point_contacts(), impulse_status.max_point_contacts());
  EXPECT_EQ(contact_status.num_active_contacts(), impulse_status.num_active_impulse());
  EXPECT_FALSE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(impulse_status.dimp(), 0);
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
}


TEST_F(DiscreteEventTest, comparison) {
  DiscreteEvent impulse_status1(max_point_contacts);
  DiscreteEvent impulse_status2(max_point_contacts);
  impulse_status1.activateImpulse({5, 6, 7});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status2.activateImpulse({1, 2, 3});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status1.activateImpulse({1, 2, 3});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status2.activateImpulse({5, 6, 7});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
}


TEST_F(DiscreteEventTest, activate) {
  DiscreteEvent impulse_status(max_point_contacts);
  ContactStatus contact_status(max_point_contacts);
  contact_status.activateContact(3);
  impulse_status.activateImpulse(3);
  EXPECT_EQ(contact_status.num_active_contacts(), impulse_status.num_active_impulse());
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3);
  EXPECT_EQ(contact_status.dimf(), impulse_status.dimp());
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_TRUE(impulse_status.isImpulseActive(3));
  for (int i=4; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  contact_status.activateContacts({5, 6});
  impulse_status.activateImpulse({5, 6});
  EXPECT_EQ(contact_status.num_active_contacts(), 3);
  EXPECT_EQ(contact_status.num_active_contacts(), impulse_status.num_active_impulse());
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 9);
  EXPECT_EQ(contact_status.dimf(), impulse_status.dimp());
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_TRUE(impulse_status.isImpulseActive(3));
  EXPECT_FALSE(contact_status.isContactActive(4));
  EXPECT_FALSE(impulse_status.isImpulseActive(4));
  EXPECT_TRUE(contact_status.isContactActive(5));
  EXPECT_TRUE(impulse_status.isImpulseActive(5));
  EXPECT_TRUE(contact_status.isContactActive(6));
  EXPECT_TRUE(impulse_status.isImpulseActive(6));
  for (int i=7; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
}


TEST_F(DiscreteEventTest, deactivate) {
  DiscreteEvent impulse_status(max_point_contacts);
  ContactStatus contact_status(max_point_contacts);
  contact_status.setContactStatus(std::vector<bool>(max_point_contacts, true));
  impulse_status.setDiscreteEvent(std::vector<bool>(max_point_contacts, true));
  EXPECT_EQ(contact_status.max_point_contacts(), max_point_contacts);
  EXPECT_EQ(contact_status.num_active_contacts(), max_point_contacts);
  EXPECT_EQ(impulse_status.max_point_contacts(), max_point_contacts);
  EXPECT_EQ(impulse_status.num_active_impulse(), max_point_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_point_contacts);
  EXPECT_EQ(impulse_status.dimp(), 3*max_point_contacts);
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  contact_status.deactivateContact(3);
  impulse_status.deactivateImpulse(3);
  EXPECT_EQ(contact_status.num_active_contacts(), max_point_contacts-1);
  EXPECT_EQ(impulse_status.num_active_impulse(), max_point_contacts-1);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_point_contacts-3);
  EXPECT_EQ(impulse_status.dimp(), 3*max_point_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_FALSE(impulse_status.isImpulseActive(3));
  for (int i=4; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  contact_status.deactivateContacts({5, 6});
  impulse_status.deactivateImpulse({5, 6});
  EXPECT_EQ(contact_status.num_active_contacts(), max_point_contacts-3);
  EXPECT_EQ(impulse_status.num_active_impulse(), max_point_contacts-3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_point_contacts-9);
  EXPECT_EQ(impulse_status.dimp(), 3*max_point_contacts-9);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_FALSE(impulse_status.isImpulseActive(3));
  EXPECT_TRUE(contact_status.isContactActive(4));
  EXPECT_TRUE(impulse_status.isImpulseActive(4));
  EXPECT_FALSE(contact_status.isContactActive(5));
  EXPECT_FALSE(impulse_status.isImpulseActive(5));
  EXPECT_FALSE(contact_status.isContactActive(6));
  EXPECT_FALSE(impulse_status.isImpulseActive(6));
  for (int i=7; i<contact_status.max_point_contacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}