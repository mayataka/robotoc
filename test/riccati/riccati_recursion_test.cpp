#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"

#include "test_helper.hpp"
#include "robot_factory.hpp"
#include "contact_sequence_factory.hpp"
#include "solution_factory.hpp"
#include "kkt_factory.hpp"
#include "cost_factory.hpp"
#include "constraints_factory.hpp"


namespace robotoc {

class RiccatiRecursionTest : public ::testing::TestWithParam<Robot> {
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

  std::shared_ptr<ContactSequence> createContactSequence(const Robot& robot) const;
  void test_riccatiRecursion(const Robot& robot) const;
  void test_computeDirection(const Robot& robot) const;

  int N, max_num_impulse, nthreads;
  double T, t, dt;
};


std::shared_ptr<ContactSequence > RiccatiRecursionTest::createContactSequence(const Robot& robot) const {
  return testhelper::CreateContactSequence(robot, N, max_num_impulse, t, 3*dt);
}


TEST_P(RiccatiRecursionTest, riccatiRecursion) {
  const auto robot = GetParam();
  auto cost = testhelper::CreateCost(robot);
  auto constraints = testhelper::CreateConstraints(robot);
  const auto contact_sequence = createContactSequence(robot);
  // DirectMultipleShooting dms(nthreads);
  aligned_vector<Robot> robots(nthreads, robot);
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, RiccatiRecursionTest, 
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