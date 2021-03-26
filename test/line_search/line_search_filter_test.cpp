#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/line_search/line_search_filter.hpp"


namespace idocp {

class LineSearchFilterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }
};


TEST_F(LineSearchFilterTest, test) {
  LineSearchFilter filter;
  EXPECT_TRUE(filter.isEmpty());
  const double cost = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double violation = std::abs(Eigen::VectorXd::Random(1)[0]);
  filter.augment(cost, violation);
  EXPECT_FALSE(filter.isEmpty());
  const double cost_accpectable = 0.5 * cost;
  const double violation_accpectable = 0.5 * violation;
  const double cost_unaccpectable = 2 * cost;
  const double violation_unaccpectable = 2 * violation;
  EXPECT_TRUE(filter.isAccepted(cost_accpectable, violation_accpectable));
  EXPECT_TRUE(filter.isAccepted(cost_unaccpectable, violation_accpectable));
  EXPECT_TRUE(filter.isAccepted(cost_accpectable, violation_unaccpectable));
  EXPECT_FALSE(filter.isAccepted(cost_unaccpectable, violation_unaccpectable));
  filter.clear();
  EXPECT_TRUE(filter.isEmpty());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}