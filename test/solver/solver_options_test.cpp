#include <vector>

#include <gtest/gtest.h>

#include "robotoc/solver/solver_options.hpp"


namespace robotoc {

class SolverOptionsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};


TEST_F(SolverOptionsTest, cout) {
  auto options = SolverOptions();
  EXPECT_NO_THROW(
    std::cout << options << std::endl;
  );
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}