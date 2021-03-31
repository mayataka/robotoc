#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"

#include "robot_factory.hpp"


namespace idocp {

class StateConstraintJacobianTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    max_num_impulse = 5;
  }

  virtual void TearDown() {
  }

  void test(const Robot& robot) const;

  int max_num_impulse;
};


void StateConstraintJacobianTest::test(const Robot& robot) const {
  StateConstraintJacobian jac(robot, max_num_impulse);
  SplitStateConstraintJacobian jac_ref(robot);
  for (int i=0; i<max_num_impulse; ++i) {
    auto impulse_status = robot.createImpulseStatus();
    impulse_status.setRandom();
    jac_ref.setImpulseStatus(impulse_status);
    jac[i].setImpulseStatus(impulse_status);
    EXPECT_TRUE(jac_ref.isApprox(jac[i]));
  }
}


TEST_F(StateConstraintJacobianTest, fixedBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test(robot);
}


TEST_F(StateConstraintJacobianTest, floatingBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}