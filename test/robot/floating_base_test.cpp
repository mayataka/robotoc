#include <string>

#include <gtest/gtest.h>
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/parsers/urdf.hpp"

#include "idocp/robot/floating_base.hpp"


namespace idocp {

class FloatingBaseTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(fixed_base_urdf, fixed_base_robot);
    pinocchio::urdf::buildModel(floating_base_urdf, floating_base_robot);
  }

  virtual void TearDown() {
  }

  std::string fixed_base_urdf, floating_base_urdf;
  pinocchio::Model fixed_base_robot, floating_base_robot;
};


TEST_F(FloatingBaseTest, fixedBaseRobot) {
  FloatingBase floating_base(fixed_base_robot);
  EXPECT_EQ(floating_base.dim_passive(), 0);
  EXPECT_TRUE(floating_base.passive_joint_indices().empty());
  EXPECT_FALSE(floating_base.hasFloatingBase());
}


TEST_F(FloatingBaseTest, floatingBaseRobot) {
  FloatingBase floating_base(floating_base_robot);
  EXPECT_EQ(floating_base.dim_passive(), 6);
  EXPECT_FALSE(floating_base.passive_joint_indices().empty());
  for (int i=0; i<floating_base.dim_passive(); ++i) {
    EXPECT_EQ(floating_base.passive_joint_indices()[i], i);
  }
  EXPECT_TRUE(floating_base.hasFloatingBase());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}