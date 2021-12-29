#include <vector>

#include <gtest/gtest.h>

#include "robotoc/robot/contact_status.hpp"


namespace robotoc {

class ContactStatusTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    max_num_contacts = 10;
  }

  virtual void TearDown() {
  }

  int max_num_contacts;
};


TEST_F(ContactStatusTest, constructor) {
  ContactStatus contact_status(max_num_contacts);
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_FALSE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 0);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  EXPECT_EQ(contact_status.contactPositions().size(), max_num_contacts);
  EXPECT_EQ(contact_status.contactSurfacesRotations().size(), max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.contactPosition(i).isZero());
    EXPECT_TRUE(contact_status.contactSurfaceRotation(i).isIdentity());
  }
}


TEST_F(ContactStatusTest, comparison) {
  ContactStatus contact_status1(max_num_contacts);
  ContactStatus contact_status2(max_num_contacts);
  contact_status1.activateContacts({5, 6, 7});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({1, 2, 3});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status1.activateContacts({1, 2, 3});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({5, 6, 7});
  EXPECT_TRUE(contact_status1 == contact_status2);
}


TEST_F(ContactStatusTest, activate) {
  ContactStatus contact_status(max_num_contacts);
  contact_status.activateContact(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3);
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  contact_status.activateContacts({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 9);
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_FALSE(contact_status.isContactActive(4));
  EXPECT_TRUE(contact_status.isContactActive(5));
  EXPECT_TRUE(contact_status.isContactActive(6));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  EXPECT_NO_THROW(
    std::cout << contact_status << std::endl;
  );
}


TEST_F(ContactStatusTest, activateAll) {
  ContactStatus contact_status(max_num_contacts);
  contact_status.activateContacts();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
}


TEST_F(ContactStatusTest, deactivate) {
  ContactStatus contact_status(max_num_contacts);
  contact_status.setActivity(std::vector<bool>(max_num_contacts, true));
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  contact_status.deactivateContact(3);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  contact_status.deactivateContacts({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-9);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  EXPECT_TRUE(contact_status.isContactActive(4));
  EXPECT_FALSE(contact_status.isContactActive(5));
  EXPECT_FALSE(contact_status.isContactActive(6));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
}


TEST_F(ContactStatusTest, deactivateAll) {
  ContactStatus contact_status(max_num_contacts);
  contact_status.setRandom();
  contact_status.deactivateContacts();
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}