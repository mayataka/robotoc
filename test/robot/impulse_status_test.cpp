#include <vector>

#include <gtest/gtest.h>

#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

class ImpulseStatusTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    max_num_contacts = 10;
    contact_types = std::vector<ContactType>(max_num_contacts, ContactType::PointContact);
    contact_frame_names.clear();
    for (int i=0; i<max_num_contacts; ++i) {
      contact_frame_names.push_back(std::to_string(i));
    }
  }

  virtual void TearDown() {
  }

  int max_num_contacts;
  std::vector<ContactType> contact_types;
  std::vector<std::string> contact_frame_names;
};


TEST_F(ImpulseStatusTest, constructor) {
  ImpulseStatus impulse_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  EXPECT_EQ(contact_status.maxNumContacts(), impulse_status.maxNumContacts());
  EXPECT_FALSE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(impulse_status.dimf(), 0);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  EXPECT_EQ(impulse_status.contactPlacements().size(), max_num_contacts);
  EXPECT_EQ(impulse_status.contactPositions().size(), max_num_contacts);
  EXPECT_EQ(impulse_status.contactRotations().size(), max_num_contacts);
  for (int i=0; i<impulse_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(impulse_status.contactPlacement(i).isIdentity());
    EXPECT_TRUE(impulse_status.contactPosition(i).isZero());
    EXPECT_TRUE(impulse_status.contactRotation(i).isIdentity());
  }
}


TEST_F(ImpulseStatusTest, comparison) {
  ImpulseStatus impulse_status1(contact_types, contact_frame_names);
  ImpulseStatus impulse_status2(contact_types, contact_frame_names);
  impulse_status1.activateImpulses({5, 6, 7});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status2.activateImpulses({1, 2, 3});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status1.activateImpulses({1, 2, 3});
  EXPECT_FALSE(impulse_status1 == impulse_status2);
  impulse_status2.activateImpulses({5, 6, 7});
  EXPECT_TRUE(impulse_status1 == impulse_status2);
}


TEST_F(ImpulseStatusTest, activate) {
  ImpulseStatus impulse_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  contact_status.activateContact(3);
  impulse_status.activateImpulse(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3);
  EXPECT_EQ(contact_status.dimf(), impulse_status.dimf());
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_TRUE(impulse_status.isImpulseActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  contact_status.activateContacts({5, 6});
  impulse_status.activateImpulses({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 9);
  EXPECT_EQ(contact_status.dimf(), impulse_status.dimf());
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
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impulse_status.isImpulseActive(i));
  }
  EXPECT_NO_THROW(
    std::cout << impulse_status << std::endl;
  );
}


TEST_F(ImpulseStatusTest, deactivate) {
  ImpulseStatus impulse_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  for (int i=0; i<contact_types.size(); ++i) {
    contact_status.activateContact(i);
    impulse_status.activateImpulse(i);
  }
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_EQ(impulse_status.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts);
  EXPECT_EQ(impulse_status.dimf(), 3*max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  contact_status.deactivateContact(3);
  impulse_status.deactivateImpulse(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-3);
  EXPECT_EQ(impulse_status.dimf(), 3*max_num_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_FALSE(impulse_status.isImpulseActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
  contact_status.deactivateContacts({5, 6});
  impulse_status.deactivateImpulses({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impulse_status.hasActiveImpulse());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-9);
  EXPECT_EQ(impulse_status.dimf(), 3*max_num_contacts-9);
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
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impulse_status.isImpulseActive(i));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}