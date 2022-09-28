#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/sto/periodic_switching_time_cost.hpp"

#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"


namespace robotoc {

class STOCostFunctionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impulse = 5;
    T = 1.0;
    t = 0.1; 

    auto robot = testhelper::CreateQuadrupedalRobot();
    kkt_matrix = KKTMatrix(robot, N, max_num_impulse);
    kkt_residual = KKTResidual(robot, N, max_num_impulse);

    const double dt = T / N;
    contact_sequence 
        = testhelper::CreateContactSequence(robot, N, max_num_impulse, 0, 3*dt);

    time_discretization = TimeDiscretization(T, N, 2*max_num_impulse);
    time_discretization.discretize(contact_sequence, t);
  }

  virtual void TearDown() {
  }

  int N, max_num_impulse;
  double T, t;
  KKTMatrix kkt_matrix;
  KKTResidual kkt_residual;
  std::shared_ptr<ContactSequence> contact_sequence;
  TimeDiscretization time_discretization;
};


TEST_F(STOCostFunctionTest, test) {
  auto cost = std::make_shared<STOCostFunction>();
  const double t_start = 0.5;
  const double period = 0.1;
  auto period_cost = std::make_shared<PeriodicSwitchingTimeCost>(period, t_start);
  cost->push_back(period_cost);
  const double cost_value1 = cost->evalCost(time_discretization);
  const double cost_value2 = cost->linearizeCost(time_discretization, kkt_residual);
  const double cost_value3 = cost->quadratizeCost(time_discretization, kkt_matrix, kkt_residual);
  EXPECT_DOUBLE_EQ(cost_value1, cost_value2);
  EXPECT_DOUBLE_EQ(cost_value2, cost_value3);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
