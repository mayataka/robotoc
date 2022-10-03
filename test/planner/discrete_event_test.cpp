#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/planner/discrete_event.hpp"


namespace robotoc {

class DiscreteEventTest : public ::testing::Test {
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
    contact_status = ContactStatus(contact_types, contact_frame_names);
  }

  virtual void TearDown() {
  }

  int max_num_contacts;
  std::vector<ContactType> contact_types;
  std::vector<std::string> contact_frame_names;
  ContactStatus contact_status;
};


TEST_F(DiscreteEventTest, constructor) {
  ContactStatus cs_before(contact_types, contact_frame_names), cs_after(contact_types, contact_frame_names);
  std::vector<Eigen::Vector3d> contact_positions;
  for (int i=0; i<max_num_contacts; ++i) {
    contact_positions.push_back(Eigen::Vector3d::Random());
  }
  cs_after.setContactPlacements(contact_positions);
  cs_before.activateContacts({1, 2, 3, 6, 7});
  cs_after.activateContacts({3, 5, 6, 7, 8, 9});
  DiscreteEvent discrete_event(cs_before, cs_after);
  EXPECT_EQ(discrete_event.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(discrete_event.preContactStatus() == cs_before);
  EXPECT_TRUE(discrete_event.postContactStatus() == cs_after);
  EXPECT_TRUE(discrete_event.existDiscreteEvent());
  EXPECT_TRUE(discrete_event.existImpulse());
  EXPECT_TRUE(discrete_event.existLift());
  for (int i=0; i<max_num_contacts; ++i) {
    EXPECT_TRUE(contact_positions[i].isApprox(discrete_event.impulseStatus().contactPositions()[i]));
  }
  EXPECT_EQ(discrete_event.eventType(), DiscreteEventType::Impact);
  EXPECT_NO_THROW(
    std::cout << discrete_event << std::endl;
  );
}


TEST_F(DiscreteEventTest, impulse) {
  ContactStatus cs_before(contact_types, contact_frame_names), cs_after(contact_types, contact_frame_names);
  cs_before.activateContacts({1, 2, 3});
  cs_after.activateContacts({1, 2, 3, 4, 5, 6});
  DiscreteEvent discrete_event(cs_before, cs_after);
  EXPECT_EQ(discrete_event.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(discrete_event.preContactStatus() == cs_before);
  EXPECT_TRUE(discrete_event.postContactStatus() == cs_after);
  EXPECT_TRUE(discrete_event.existDiscreteEvent());
  EXPECT_TRUE(discrete_event.existImpulse());
  EXPECT_FALSE(discrete_event.existLift());
  EXPECT_EQ(discrete_event.eventType(), DiscreteEventType::Impact);
  EXPECT_NO_THROW(
    std::cout << discrete_event << std::endl;
  );
}


TEST_F(DiscreteEventTest, lift) {
  ContactStatus cs_before(contact_types, contact_frame_names), cs_after(contact_types, contact_frame_names);
  cs_before.activateContacts({1, 2, 3, 4, 5, 6});
  cs_after.activateContacts({1, 2, 3});
  DiscreteEvent discrete_event(cs_before, cs_after);
  EXPECT_EQ(discrete_event.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(discrete_event.preContactStatus() == cs_before);
  EXPECT_TRUE(discrete_event.postContactStatus() == cs_after);
  EXPECT_TRUE(discrete_event.existDiscreteEvent());
  EXPECT_FALSE(discrete_event.existImpulse());
  EXPECT_TRUE(discrete_event.existLift());
  EXPECT_EQ(discrete_event.eventType(), DiscreteEventType::Lift);
  EXPECT_NO_THROW(
    std::cout << discrete_event << std::endl;
  );
}


TEST_F(DiscreteEventTest, impulseAndLift) {
  ContactStatus cs_before(contact_types, contact_frame_names), cs_after(contact_types, contact_frame_names);
  cs_before.activateContacts({1, 2, 3, 6, 7});
  cs_after.activateContacts({3, 5, 6, 7, 8, 9});
  DiscreteEvent discrete_event(cs_before, cs_after);
  EXPECT_EQ(discrete_event.maxNumContacts(), max_num_contacts);
  EXPECT_TRUE(discrete_event.preContactStatus() == cs_before);
  EXPECT_TRUE(discrete_event.postContactStatus() == cs_after);
  EXPECT_TRUE(discrete_event.existDiscreteEvent());
  EXPECT_TRUE(discrete_event.existImpulse());
  EXPECT_TRUE(discrete_event.existLift());
  EXPECT_EQ(discrete_event.eventType(), DiscreteEventType::Impact);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}