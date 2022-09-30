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
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class LineSearchTest : public ::testing::TestWithParam<LineSearchSettings> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
    step_size_reduction_rate = 0.75;
    min_step_size = 0.05;
  }

  virtual void TearDown() {
  }

  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence) const;
  Direction createDirection(const Robot& robot) const;
  Direction createDirection(const Robot& robot, 
                            const std::shared_ptr<ContactSequence>& contact_sequence) const;
  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;

  void test(const Robot& robot, const LineSearchSettings& settings) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt, step_size_reduction_rate, min_step_size;
};


Solution LineSearchTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution LineSearchTest::createSolution(const Robot& robot, 
                                        const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


Direction LineSearchTest::createDirection(const Robot& robot) const {
  return testhelper::CreateDirection(robot, N, max_num_impulse);
}


Direction LineSearchTest::createDirection(const Robot& robot, 
                                          const std::shared_ptr<ContactSequence>& contact_sequence) const {
  return testhelper::CreateDirection(robot, contact_sequence, T, N, max_num_impulse, t);
}


std::shared_ptr<ContactSequence> LineSearchTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


void LineSearchTest::test(const Robot& robot, const LineSearchSettings& settings) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence);
  const auto d = createDirection(robot, contact_sequence);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots(nthreads, robot);
  // auto ocp = OCP(robot, cost, constraints, contact_sequence, T, N);
  OCPDef ocp;
  ocp.robot = robot;
  ocp.cost = cost;
  ocp.constraints = constraints;
  ocp.N = N;
  ocp.T = T;
  // ocp.discretize(t);
  DirectMultipleShooting dms(ocp, nthreads);
  dms.initConstraints(ocp, robots, contact_sequence, s);
  LineSearch line_search(ocp, nthreads, settings);
  EXPECT_TRUE(line_search.isFilterEmpty());
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(ocp, robots, contact_sequence, q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  if (settings.line_search_method == LineSearchMethod::Filter) {
    EXPECT_FALSE(line_search.isFilterEmpty());
  }
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_DOUBLE_EQ(line_search.computeStepSize(ocp, robots, contact_sequence, q, v, s, d, very_small_max_primal_step_size),
                   min_step_size);
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