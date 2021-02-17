#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"


namespace idocp {

class StateConstraintJacobianTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    max_num_impulse = 5;
  }

  virtual void TearDown() {
  }


  void test(const Robot& robot) const;

  std::string fixed_base_urdf, floating_base_urdf;
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
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  test(robot);
}


TEST_F(StateConstraintJacobianTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}