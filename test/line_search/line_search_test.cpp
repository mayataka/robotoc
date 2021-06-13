#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/direct_multiple_shooting.hpp"
#include "idocp/line_search/line_search_filter.hpp"
#include "idocp/line_search/line_search.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "direction_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace idocp {

class LineSearchTest : public ::testing::Test {
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
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence) const;
  Direction createDirection(const Robot& robot) const;
  Direction createDirection(const Robot& robot, const ContactSequence& contact_sequence) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void test(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt, step_size_reduction_rate, min_step_size;
};


Solution LineSearchTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution LineSearchTest::createSolution(const Robot& robot, 
                                        const ContactSequence& contact_sequence) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
}


Direction LineSearchTest::createDirection(const Robot& robot) const {
  return testhelper::CreateDirection(robot, N, max_num_impulse);
}


Direction LineSearchTest::createDirection(const Robot& robot, 
                                          const ContactSequence& contact_sequence) const {
  return testhelper::CreateDirection(robot, contact_sequence, T, N, max_num_impulse, t);
}


ContactSequence LineSearchTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


void LineSearchTest::test(const Robot& robot) const {
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
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  DirectMultipleShooting dms(N, max_num_impulse, nthreads);
  dms.initConstraints(ocp, robots, contact_sequence, s);
  LineSearch line_search(robot, N, max_num_impulse, nthreads);
  EXPECT_TRUE(line_search.isFilterEmpty());
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(ocp, robots, contact_sequence, q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  EXPECT_FALSE(line_search.isFilterEmpty());
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_DOUBLE_EQ(line_search.computeStepSize(ocp, robots, contact_sequence, q, v, s, d, very_small_max_primal_step_size),
                   min_step_size);
}


TEST_F(LineSearchTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot();
  test(robot);
  robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(LineSearchTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot();
  test(robot);
  robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}