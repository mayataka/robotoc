#include <gtest/gtest.h>

#include "idocp/constraints/constraints_data.hpp"

namespace idocp {

class ConstraintsDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }

};


TEST_F(ConstraintsDataTest, defaultConstructor) {
  ConstraintsData data;
  EXPECT_FALSE(data.isPositionLevelValid());
  EXPECT_FALSE(data.isVelocityLevelValid());
  EXPECT_FALSE(data.isAccelerationLevelValid());
  EXPECT_FALSE(data.isImpulseLevelValid());
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_TRUE(data.acceleration_level_data.empty());
  EXPECT_TRUE(data.impulse_level_data.empty());
}


TEST_F(ConstraintsDataTest, timestep0) {
  const int time_step = 0;
  ConstraintsData data(time_step);
  EXPECT_FALSE(data.isPositionLevelValid());
  EXPECT_FALSE(data.isVelocityLevelValid());
  EXPECT_TRUE(data.isAccelerationLevelValid());
  EXPECT_FALSE(data.isImpulseLevelValid());
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_TRUE(data.acceleration_level_data.empty());
  EXPECT_TRUE(data.impulse_level_data.empty());
}


TEST_F(ConstraintsDataTest, timestep1) {
  const int time_step = 1;
  ConstraintsData data(time_step);
  EXPECT_FALSE(data.isPositionLevelValid());
  EXPECT_TRUE(data.isVelocityLevelValid());
  EXPECT_TRUE(data.isAccelerationLevelValid());
  EXPECT_FALSE(data.isImpulseLevelValid());
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_TRUE(data.acceleration_level_data.empty());
  EXPECT_TRUE(data.impulse_level_data.empty());
}


TEST_F(ConstraintsDataTest, timestep2) {
  const int time_step = 2;
  ConstraintsData data(time_step);
  EXPECT_TRUE(data.isPositionLevelValid());
  EXPECT_TRUE(data.isVelocityLevelValid());
  EXPECT_TRUE(data.isAccelerationLevelValid());
  EXPECT_FALSE(data.isImpulseLevelValid());
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_TRUE(data.acceleration_level_data.empty());
  EXPECT_TRUE(data.impulse_level_data.empty());
}


TEST_F(ConstraintsDataTest, timestepminus1) {
  const int time_step = -1;
  ConstraintsData data(time_step);
  EXPECT_FALSE(data.isPositionLevelValid());
  EXPECT_FALSE(data.isVelocityLevelValid());
  EXPECT_FALSE(data.isAccelerationLevelValid());
  EXPECT_TRUE(data.isImpulseLevelValid());
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_TRUE(data.acceleration_level_data.empty());
  EXPECT_TRUE(data.impulse_level_data.empty());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}