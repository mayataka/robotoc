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
    contact_types = std::vector<ContactType>(max_num_contacts, ContactType::PointContact);
    contact_frame_names.clear();
    for (int i=0; i<max_num_contacts; ++i) {
      contact_frame_names.push_back(std::to_string(i));
    }
  }

  virtual void TearDown() {
  }

  void checkContactStatusAvtiveByIndexAndByName(const ContactStatus& contact_status) const {
    for (int i=0; i<contact_status.maxNumContacts(); ++i) {
      EXPECT_EQ(contact_status.isContactActive(i), 
                contact_status.isContactActive(contact_frame_names[i]));
    }
  }

  int max_num_contacts;
  std::vector<ContactType> contact_types;
  std::vector<std::string> contact_frame_names;
};


TEST_F(ContactStatusTest, constructor) {
  ContactStatus contact_status(contact_types, contact_frame_names);
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_FALSE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 0);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
    EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[i]));
  }
  EXPECT_EQ(contact_status.contactPlacements().size(), max_num_contacts);
  EXPECT_EQ(contact_status.contactPositions().size(), max_num_contacts);
  EXPECT_EQ(contact_status.contactRotations().size(), max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.contactPlacement(i).isIdentity());
    EXPECT_TRUE(contact_status.contactPosition(i).isZero());
    EXPECT_TRUE(contact_status.contactRotation(i).isIdentity());
  }
  EXPECT_EQ(contact_status.frictionCoefficients().size(), max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_DOUBLE_EQ(contact_status.frictionCoefficient(i), 0.7); // default value
  }
}


TEST_F(ContactStatusTest, comparison) {
  ContactStatus contact_status1(contact_types, contact_frame_names);
  ContactStatus contact_status2(contact_types, contact_frame_names);
  contact_status1.activateContacts({5, 6, 7});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({1, 2, 3});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status1.activateContacts({1, 2, 3});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({5, 6, 7});
  EXPECT_TRUE(contact_status1 == contact_status2);
}


TEST_F(ContactStatusTest, comparison2) {
  ContactStatus contact_status1(contact_types, contact_frame_names);
  ContactStatus contact_status2(contact_types, contact_frame_names);
  contact_status1.activateContacts({"5", "6", "7"});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({"1", "2", "3"});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status1.activateContacts({"1", "2", "3"});
  EXPECT_FALSE(contact_status1 == contact_status2);
  contact_status2.activateContacts({"5", "6", "7"});
  EXPECT_TRUE(contact_status1 == contact_status2);
}


TEST_F(ContactStatusTest, activate) {
  ContactStatus contact_status(contact_types, contact_frame_names);
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
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  contact_status.activateContacts({5, 6});
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 9);
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_TRUE(contact_status.isContactActive(3));
  EXPECT_FALSE(contact_status.isContactActive(4));
  EXPECT_TRUE(contact_status.isContactActive(5));
  EXPECT_TRUE(contact_status.isContactActive(6));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(i));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_NO_THROW(
    std::cout << contact_status << std::endl;
  );
}


TEST_F(ContactStatusTest, activate2) {
  ContactStatus contact_status(contact_types, contact_frame_names);
  contact_status.activateContact(contact_frame_names[3]);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3);
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[i]));
  }
  EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[3]));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[i]));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  contact_status.activateContacts(std::vector<std::string>({"5", "6"}));
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 9);
  for (int i=0; i<3; ++i) {
    EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[i]));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[3]));
  EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[4]));
  EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[5]));
  EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[6]));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[i]));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_NO_THROW(
    std::cout << contact_status << std::endl;
  );
}


TEST_F(ContactStatusTest, deactivate) {
  ContactStatus contact_status(contact_types, contact_frame_names);
  for (int i=0; i<contact_types.size(); ++i) {
    contact_status.activateContact(i);
  }
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  contact_status.deactivateContact(3);
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  EXPECT_FALSE(contact_status.isContactActive(3));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(i));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
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
  checkContactStatusAvtiveByIndexAndByName(contact_status);
}


TEST_F(ContactStatusTest, deactivate2) {
  ContactStatus contact_status(contact_types, contact_frame_names);
  for (int i=0; i<contact_types.size(); ++i) {
    contact_status.activateContact(contact_frame_names[i]);
  }
  EXPECT_EQ(contact_status.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts);
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[i]));
  }
  contact_status.deactivateContact(contact_frame_names[3]);
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-3);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[i]));
  }
  EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[3]));
  for (int i=4; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[i]));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
  contact_status.deactivateContacts(std::vector<std::string>({"5", "6"}));
  EXPECT_TRUE(contact_status.hasActiveContacts());
  EXPECT_EQ(contact_status.dimf(), 3*max_num_contacts-9);
  for (int i=0; i<3; ++i) {
    EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[i]));
  }
  EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[3]));
  EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[4]));
  EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[5]));
  EXPECT_FALSE(contact_status.isContactActive(contact_frame_names[6]));
  for (int i=7; i<contact_status.maxNumContacts(); ++i) {
    EXPECT_TRUE(contact_status.isContactActive(contact_frame_names[i]));
  }
  checkContactStatusAvtiveByIndexAndByName(contact_status);
}
} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}