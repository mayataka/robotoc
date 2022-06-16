#include <vector>

#include <gtest/gtest.h>

#include "robotoc/mpc/moving_window_filter.hpp"


namespace robotoc {

class MovingWindowFilterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(MovingWindowFilterTest, test) {
  const double time_length = 0.5;
  const double min_sampling_period = 0.1;
  MovingWindowFilter<2> filter(time_length, min_sampling_period);
  double t = 0;
  for (int i=0; i<100; ++i) {
    EXPECT_NO_THROW(
      filter.push_back(t, Eigen::Vector2d::Random());
    );
    t += 0.05;
  }
  EXPECT_EQ(filter.size(), 4);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}