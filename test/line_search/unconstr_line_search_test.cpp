#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
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
    time_discretization.resize(N+1);
    for (int i=0; i<=N; ++i) {
      time_discretization[i].t = t + dt * i;
      time_discretization[i].dt = dt;
      time_discretization[i].phase = -1;
      time_discretization[i].stage = i;
      time_discretization[i].impact_index = -1;
      time_discretization[i].lift_index = -1;
      time_discretization[i].stage_in_phase = i;
      time_discretization[i].num_grids_in_phase = N;
    }
  }

  virtual void TearDown() {
  }

  Robot robot;
  int dimv, N, nthreads;
  double T, dt, t, step_size_reduction_rate, min_step_size;
  std::vector<GridInfo> time_discretization;
};


TEST_F(UnconstrLineSearchTest, UnconstrOCP) {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto s = testhelper::CreateSolution(robot, N);
  const auto d = testhelper::CreateDirection(robot, N);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  aligned_vector<Robot> robots(nthreads, robot);
  OCP ocp(robot, cost, constraints, T, N);
  UnconstrDirectMultipleShooting dms(ocp, nthreads);
  dms.initConstraints(robots, time_discretization, s);
  UnconstrLineSearch line_search(ocp);
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(dms, robots, time_discretization, 
                                                       q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size2 = line_search.computeStepSize(dms, robots, time_discretization, 
                                                        q, v, s, d, max_primal_step_size);
  EXPECT_DOUBLE_EQ(step_size2, min_step_size);
}


TEST_F(UnconstrLineSearchTest, UnconstrParNMPC) {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto s = testhelper::CreateSolution(robot, N);
  const auto d = testhelper::CreateDirection(robot, N);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  aligned_vector<Robot> robots(nthreads, robot);
  OCP ocp(robot, cost, constraints, T, N);
  UnconstrBackwardCorrection backward_correction(ocp, nthreads);
  backward_correction.initConstraints(robots, time_discretization, s);
  UnconstrLineSearch line_search(ocp);
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(backward_correction, robots, time_discretization, 
                                                       q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  const double step_size2 = line_search.computeStepSize(backward_correction, robots, time_discretization, 
                                                        q, v, s, d, max_primal_step_size);
  EXPECT_DOUBLE_EQ(step_size2, min_step_size);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}