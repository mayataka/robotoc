#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/line_search/line_search_filter.hpp"
#include "robotoc/line_search/line_search_settings.hpp"
#include "robotoc/line_search/line_search.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "direction_factory.hpp"
#include "kkt_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class LineSearchTest : public ::testing::TestWithParam<LineSearchSettings> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impact = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
    step_size_reduction_rate = 0.75;
  }

  virtual void TearDown() {
  }
  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;

  void test(const Robot& robot, const LineSearchSettings& settings) const;

  int N, max_num_impact, nthreads;
  double T, t, dt, step_size_reduction_rate;
};


std::shared_ptr<ContactSequence> LineSearchTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impact, t, 3*dt);
}


void LineSearchTest::test(const Robot& robot, const LineSearchSettings& settings) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  TimeDiscretization time_discretization(T, N);
  time_discretization.discretize(contact_sequence, t);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_residual = testhelper::CreateKKTResidual(robot, contact_sequence, time_discretization);
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots(nthreads, robot);
  OCP ocp;
  ocp.robot = robot;
  ocp.cost = cost;
  ocp.constraints = constraints;
  ocp.contact_sequence = contact_sequence;
  ocp.N = N;
  ocp.T = T;
  DirectMultipleShooting dms(ocp, nthreads);
  const auto s = testhelper::CreateSolution(robot, contact_sequence, time_discretization);
  const auto d = testhelper::CreateDirection(robot, contact_sequence, time_discretization);
  dms.initConstraints(robots, time_discretization, s);
  dms.evalOCP(robots, time_discretization, q, v, s, kkt_residual);
  LineSearch line_search(ocp, settings);
  const double max_primal_step_size = settings.min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-settings.min_step_size);
  const double step_size = line_search.computeStepSize(dms, robots, time_discretization, 
                                                       q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size > 0.0);
  const double very_small_max_primal_step_size = settings.min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_DOUBLE_EQ(line_search.computeStepSize(dms, robots, time_discretization, q, v, s, d, very_small_max_primal_step_size),
                   very_small_max_primal_step_size);
  EXPECT_NO_THROW(
    std::cout << settings << std::endl;
  );
}


TEST_P(LineSearchTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator();
  auto settings = GetParam();
  test(robot, settings);
  robot = testhelper::CreateRobotManipulator(dt);
  test(robot, settings);
}


TEST_P(LineSearchTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot();
  auto settings = GetParam();
  test(robot, settings);
  robot = testhelper::CreateQuadrupedalRobot(dt);
  test(robot, settings);
}


auto createFilterSettings = []() {
  LineSearchSettings filter_settings;
  filter_settings.line_search_method = LineSearchMethod::Filter;
  return filter_settings;
};

auto createBacktrackSettings = []() {
  LineSearchSettings backtrack_settings;
  backtrack_settings.line_search_method = LineSearchMethod::MeritBacktracking;
  return backtrack_settings;
};

INSTANTIATE_TEST_SUITE_P(ParamtererizedTest, LineSearchTest, 
                         ::testing::Values(createFilterSettings(), createBacktrackSettings()));

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}