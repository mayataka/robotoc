#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/sto/sto_cost_function.hpp"

#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"


namespace robotoc {

class STOCostFunctionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impact = 5;
    T = 1.0;
    t = 0.1; 

    auto robot = testhelper::CreateQuadrupedalRobot();

    const double dt = T / N;
    contact_sequence 
        = testhelper::CreateContactSequence(robot, N, max_num_impact, 0, 3*dt);

    time_discretization = TimeDiscretization(T, N, 2*max_num_impact);
    time_discretization.discretize(contact_sequence, t);
  }

  virtual void TearDown() {
  }

  int N, max_num_impact;
  double T, t;
  Eigen::VectorXd lt;
  Eigen::MatrixXd Qtt;
  std::shared_ptr<ContactSequence> contact_sequence;
  TimeDiscretization time_discretization;
};


TEST_F(STOCostFunctionTest, test) {
  auto cost = std::make_shared<STOCostFunction>();

  const double cost_value1 = cost->evalCost(time_discretization);
  const double cost_value2 = cost->linearizeCost(time_discretization, lt);
  const double cost_value3 = cost->quadratizeCost(time_discretization, lt, Qtt);
  EXPECT_DOUBLE_EQ(cost_value1, cost_value2);
  EXPECT_DOUBLE_EQ(cost_value2, cost_value3);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
