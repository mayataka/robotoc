#include <vector>

#include <gtest/gtest.h>

#include "robotoc/mpc/raibert_heuristic.hpp"


namespace robotoc {

class RaibertHeuristicTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(RaibertHeuristicTest, test) {
  const double stance_time = 0.5;
  const double gain = 0.2;
  RaibertHeuristic heuristic(stance_time, gain);
  Eigen::Vector2d vcom; 
  Eigen::Vector2d vcom_cmd;
  vcom << 0.1, 0.0;
  vcom_cmd << 0.2, 0.0;
  EXPECT_NO_THROW(
    heuristic.planStepLength(vcom, vcom_cmd, 0.0);
  );
  std::cout << heuristic.stepLength().transpose() << std::endl;
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}