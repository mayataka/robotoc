#include <vector>

#include <gtest/gtest.h>

#include "robotoc/solver/solver_statistics.hpp"


namespace robotoc {

class SolverStatisticsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(SolverStatisticsTest, cout) {
  SolverStatistics statistics;
  EXPECT_NO_THROW(
    std::cout << statistics << std::endl;
  );
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}