#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class DirectMultipleShootingTest : public ::testing::TestWithParam<Robot> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    N = 20;
    max_num_impulse = 5;
    nthreads = 4;
    T = 1;
    t = std::abs(Eigen::VectorXd::Random(1)[0]);
    dt = T / N;
  }

  virtual void TearDown() {
  }

  // Solution createSolution(const Robot& robot) const;
  // Solution createSolution(const Robot& robot, 
  //                         const std::shared_ptr<ContactSequence>& contact_sequence) const;
  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


// Solution DirectMultipleShootingTest::createSolution(const Robot& robot) const {
//   return testhelper::CreateSolution(robot, N, max_num_impulse);
// }


// Solution DirectMultipleShootingTest::createSolution(const Robot& robot, 
//                                                     const std::shared_ptr<ContactSequence>& contact_sequence) const {
//   return testhelper::CreateSolution(robot, contact_sequence, T, N, max_num_impulse, t);
// }


std::shared_ptr<ContactSequence> DirectMultipleShootingTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


TEST_P(DirectMultipleShootingTest, evalKKT) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  TimeDiscretization time_discretization(T, N, 3*max_num_impulse);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const auto s = testhelper::CreateSolution(robot, contact_sequence, time_discretization);
  aligned_vector<Robot> robots(nthreads, robot);
  // auto kkt_matrix = testhelper::C
  // auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  // auto kkt_matrix_ref = kkt_matrix;
  // auto kkt_residual_ref = kkt_residual;
}


TEST_P(DirectMultipleShootingTest, integrateSolution) {
  auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  TimeDiscretization time_discretization(T, N, 3*max_num_impulse);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  const Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  const auto s = testhelper::CreateSolution(robot, contact_sequence, time_discretization);
  aligned_vector<Robot> robots(nthreads, robot);
  // auto kkt_matrix = testhelper::C
  // auto kkt_residual = KKTResidual(robot, N, max_num_impulse);
  // auto kkt_matrix_ref = kkt_matrix;
  // auto kkt_residual_ref = kkt_residual;
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, DirectMultipleShootingTest, 
  ::testing::Values(testhelper::CreateRobotManipulator(),
                    testhelper::CreateRobotManipulator(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateQuadrupedalRobot(),
                    testhelper::CreateQuadrupedalRobot(std::abs(Eigen::VectorXd::Random(1)[0])),
                    testhelper::CreateHumanoidRobot(),
                    testhelper::CreateHumanoidRobot(std::abs(Eigen::VectorXd::Random(1)[0])))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
