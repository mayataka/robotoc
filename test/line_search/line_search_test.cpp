#include <string>
#include <memory>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/ocp_linearizer.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"
#include "idocp/line_search/line_search_filter.hpp"
#include "idocp/line_search/line_search.hpp"

#include "test_helper.hpp"


namespace idocp {

class LineSearchTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    N = 20;
    max_num_impulse = 5;
    nthreads = 1;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dtau = T / N;
    step_size_reduction_rate = 0.75;
    min_step_size = 0.05;
  }

  virtual void TearDown() {
  }

  Solution createSolution(const Robot& robot) const;
  Solution createSolution(const Robot& robot, const ContactSequence& contact_sequence, const bool is_parnmpc) const;
  Direction createDirection(const Robot& robot) const;
  Direction createDirection(const Robot& robot, const ContactSequence& contact_sequence, const bool is_parnmpc) const;
  ContactSequence createContactSequence(const Robot& robot) const;

  void testOCP(const Robot& robot) const;
  void testParNMPC(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
  int N, max_num_impulse, nthreads;
  double T, t, dtau, step_size_reduction_rate, min_step_size;
  std::shared_ptr<CostFunction> cost;
  std::shared_ptr<Constraints> constraints;
};


Solution LineSearchTest::createSolution(const Robot& robot) const {
  return testhelper::CreateSolution(robot, N, max_num_impulse);
}


Solution LineSearchTest::createSolution(const Robot& robot, 
                                        const ContactSequence& contact_sequence, 
                                        const bool is_parnmpc) const {
  return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t, is_parnmpc);
}


Direction LineSearchTest::createDirection(const Robot& robot) const {
  return testhelper::CreateDirection(robot, N, max_num_impulse);
}


Direction LineSearchTest::createDirection(const Robot& robot, 
                                          const ContactSequence& contact_sequence,
                                          const bool is_parnmpc) const {
  return testhelper::CreateDirection(robot, contact_sequence, T, N, max_num_impulse, t, is_parnmpc);
}


ContactSequence LineSearchTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dtau);
}


void LineSearchTest::testOCP(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence, false);
  const auto d = createDirection(robot, contact_sequence, false);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot> robots(nthreads, robot);
  auto ocp = OCP(robot, cost, constraints, T, N, max_num_impulse);
  ocp.discretize(contact_sequence, t);
  OCPLinearizer linearizer(N, max_num_impulse, nthreads);
  linearizer.initConstraints(ocp, robots, contact_sequence, s);
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


void LineSearchTest::testParNMPC(const Robot& robot) const {
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  const auto s = createSolution(robot, contact_sequence, true);
  const auto d = createDirection(robot, contact_sequence, true);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  auto kkt_residual_ref = kkt_residual;
  std::vector<Robot> robots(nthreads, robot);
  auto parnmpc = ParNMPC(robot, cost, constraints, T, N, max_num_impulse);
  parnmpc.discretize(contact_sequence, t);
  ParNMPCLinearizer linearizer(N, max_num_impulse, nthreads);
  linearizer.initConstraints(parnmpc, robots, contact_sequence, s);
  LineSearch line_search(robot, N, max_num_impulse, nthreads);
  EXPECT_TRUE(line_search.isFilterEmpty());
  const double max_primal_step_size = min_step_size + std::abs(Eigen::VectorXd::Random(1)[0]) * (1-min_step_size);
  const double step_size = line_search.computeStepSize(parnmpc, robots, contact_sequence, q, v, s, d, max_primal_step_size);
  EXPECT_TRUE(step_size <= max_primal_step_size);
  EXPECT_TRUE(step_size >= min_step_size);
  EXPECT_FALSE(line_search.isFilterEmpty());
  const double very_small_max_primal_step_size = min_step_size * std::abs(Eigen::VectorXd::Random(1)[0]);
  EXPECT_DOUBLE_EQ(line_search.computeStepSize(parnmpc, robots, contact_sequence, q, v, s, d, very_small_max_primal_step_size),
                   min_step_size);
}


TEST_F(LineSearchTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testOCP(robot);
  testParNMPC(robot);
  std::vector<int> contact_frames = {18};
  robot = Robot(fixed_base_urdf, contact_frames);
  testOCP(robot);
  testParNMPC(robot);
}


TEST_F(LineSearchTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testOCP(robot);
  testParNMPC(robot);
  std::vector<int> contact_frames = {14, 24, 34, 44};
  robot = Robot(floating_base_urdf, contact_frames);
  testOCP(robot);
  testParNMPC(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}