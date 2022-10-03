#include <vector>

#include <gtest/gtest.h>

#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impact_status.hpp"


namespace robotoc {

class ImpactStatusTest : public ::testing::Test {
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


TEST_F(ImpactStatusTest, constructor) {
  ImpactStatus impact_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  EXPECT_EQ(contact_status.maxNumContacts(), impact_status.maxNumContacts());
  EXPECT_FALSE(impact_status.hasActiveImpact());
  EXPECT_EQ(impact_status.dimf(), 0);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(impact_status.isImpactActive(i));
  }
  EXPECT_EQ(impact_status.contactPlacements().size(), max_num_contacts);
  EXPECT_EQ(impact_status.contactPositions().size(), max_num_contacts);
  EXPECT_EQ(impact_status.contactRotations().size(), max_num_contacts);
  for (int i=0; i<impact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(impact_status.contactPlacement(i).isIdentity());
    EXPECT_TRUE(impact_status.contactPosition(i).isZero());
    EXPECT_TRUE(impact_status.contactRotation(i).isIdentity());
  }
}


TEST_F(ImpactStatusTest, comparison) {
  ImpactStatus impact_status1(contact_types, contact_frame_names);
  ImpactStatus impact_status2(contact_types, contact_frame_names);
  impact_status1.activateImpacts({5, 6, 7});
  EXPECT_FALSE(impact_status1 == impact_status2);
  impact_status2.activateImpacts({1, 2, 3});
  EXPECT_FALSE(impact_status1 == impact_status2);
  impact_status1.activateImpacts({1, 2, 3});
  EXPECT_FALSE(impact_status1 == impact_status2);
  impact_status2.activateImpacts({5, 6, 7});
  EXPECT_TRUE(impact_status1 == impact_status2);
}


TEST_F(ImpactStatusTest, activate) {
  ImpactStatus impact_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  contact_status.activateContact(3);
  impact_status.activateImpact(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impact_status.hasActiveImpact());
  EXPECT_EQ(contact_status.dimf(), 3);
  EXPECT_EQ(contact_status.dimf(), impact_status.dimf());
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impact_status.isImpactActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_TRUE(impact_status.isImpactActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impact_status.isImpactActive(i));
  }
  contact_status.activateContacts({5, 6});
  impact_status.activateImpacts({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impact_status.hasActiveImpact());
  EXPECT_EQ(contact_status.dimf(), 9);
  EXPECT_EQ(contact_status.dimf(), impact_status.dimf());
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impact_status.isImpactActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_TRUE(impact_status.isImpactActive(3));
  EXPECT_FALSE(contact_status.isContactActive(4));
  EXPECT_FALSE(impact_status.isImpactActive(4));
  EXPECT_TRUE(contact_status.isContactActive(5));
  EXPECT_TRUE(impact_status.isImpactActive(5));
  EXPECT_TRUE(contact_status.isContactActive(6));
  EXPECT_TRUE(impact_status.isImpactActive(6));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(impact_status.isImpactActive(i));
  }
  EXPECT_NO_THROW(
    std::cout << impact_status << std::endl;
  );
}


TEST_F(ImpactStatusTest, deactivate) {
  ImpactStatus impact_status(contact_types, contact_frame_names);
  ContactStatus contact_status(contact_types, contact_frame_names);
  for (int i=0; i<contact_types.size(); ++i) {
    contact_status.activateContact(i);
    impact_status.activateImpact(i);
  }
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_EQ(impact_status.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impact_status.hasActiveImpact());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts);
  EXPECT_EQ(impact_status.dimf(), 3*max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impact_status.isImpactActive(i));
  }
  contact_status.deactivateContact(3);
  impact_status.deactivateImpact(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impact_status.hasActiveImpact());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-3);
  EXPECT_EQ(impact_status.dimf(), 3*max_num_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impact_status.isImpactActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_FALSE(impact_status.isImpactActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impact_status.isImpactActive(i));
  }
  contact_status.deactivateContacts({5, 6});
  impact_status.deactivateImpacts({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_TRUE(impact_status.hasActiveImpact());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-9);
  EXPECT_EQ(impact_status.dimf(), 3*max_num_contacts-9);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impact_status.isImpactActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_FALSE(impact_status.isImpactActive(3));
  EXPECT_TRUE(contact_status.isContactActive(4));
  EXPECT_TRUE(impact_status.isImpactActive(4));
  EXPECT_FALSE(contact_status.isContactActive(5));
  EXPECT_FALSE(impact_status.isImpactActive(5));
  EXPECT_FALSE(contact_status.isContactActive(6));
  EXPECT_FALSE(impact_status.isImpactActive(6));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
    EXPECT_TRUE(impact_status.isImpactActive(i));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}