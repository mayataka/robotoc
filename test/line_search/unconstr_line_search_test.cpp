#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/unconstr/unconstr_ocp.hpp"
#include "robotoc/unconstr/unconstr_parnmpc.hpp"
#include "robotoc/line_search/line_search_filter.hpp"
#include "robotoc/line_search/unconstr_line_search.hpp"

#include "robot_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"
#include "solution_factory.hpp"
#include "direction_factory.hpp"


namespace robotoc {

class UnconstrLineSearchTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateRobotManipulator();
    dimv = robot.dimv();
    N = 20;
    nthreads = 4;
    T = 1;
    dt = T / N;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    step_size_reduction_rate = 0.75;
    min_step_size = 0.05;
  }

  virtual void TearDown() {
  }

  Robot robot;
  int dimv, N, nthreads;
  double T, dt, t, step_size_reduction_rate, min_step_size;
};


TEST_F(UnconstrLineSearchTest, UnconstrOCP) {
  // auto cost = testhelper::CreateCost(robot);
  // auto constraints = testhelper::CreateConstraints(robot);
  // const auto s = testhelper::CreateSolution(robot, N);
  // const auto d = testhelper::CreateDirection(robot, N);
  // const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  // const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  // aligned_vector<Robot> robots(nthreads, robot);
  // auto ocp = UnconstrOCP(robot, cost, constraints, T, N);
  // for (int i=0; i<N; ++i) {
  //   ocp[i].initConstraints(robot, i, s[i]);
  // }
  // UnconstrLineSearch line_search(ocp, nthreads);
  // EXPECT_TRUE(line_search.isFilterEmpty());
  // const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  // const double step_size = line_search.computeStepSize(ocp, robots, t, q, v, s, d, max_primal_step_size);
  // EXPECT_TRUE(step_size <= max_primal_step_size);
  // EXPECT_TRUE(step_size >= min_step_size);
  // EXPECT_FALSE(line_search.isFilterEmpty());
  // const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  // EXPECT_DOUBLE_EQ(line_search.computeStepSize(ocp, robots, t, q, v, s, d, very_small_max_primal_step_size),
  //                  min_step_size);
}


TEST_F(UnconstrLineSearchTest, UnconstrParNMPC) {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto s = testhelper::CreateSolution(robot, N);
  const auto d = testhelper::CreateDirection(robot, N);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  aligned_vector<Robot> robots(nthreads, robot);
  auto parnmpc = UnconstrParNMPC(robot, cost, constraints, T, N);
  for (int i=0; i<N-1; ++i) {
    parnmpc[i].initConstraints(robot, i+1, s[i]);
  }
  parnmpc.terminal.initConstraints(robot, N, s[N-1]);
  UnconstrLineSearch line_search(parnmpc, nthreads);
  EXPECT_TRUE(line_search.isFilterEmpty());
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(parnmpc, robots, t, q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  EXPECT_FALSE(line_search.isFilterEmpty());
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_DOUBLE_EQ(line_search.computeStepSize(parnmpc, robots, t, q, v, s, d, very_small_max_primal_step_size),
                   min_step_size);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}