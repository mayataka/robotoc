#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {

class SplitStateConstraintJacobianTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }


  static void test(const Robot& robot, const ImpulseStatus& impulse_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitStateConstraintJacobianTest::test(const Robot& robot, 
                                            const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimi = impulse_status.dimf();
  SplitStateConstraintJacobian jac(robot);
  jac.setImpulseStatus(impulse_status);
  EXPECT_EQ(jac.dimi(), dimi);
  if (robot.hasFloatingBase()) {
    EXPECT_EQ(jac.dintegrate_dq.rows(), dimv);
    EXPECT_EQ(jac.dintegrate_dq.cols(), dimv);
    EXPECT_EQ(jac.dintegrate_dv.rows(), dimv);
    EXPECT_EQ(jac.dintegrate_dv.cols(), dimv);
  }
  else {
    EXPECT_EQ(jac.dintegrate_dq.rows(), 0);
    EXPECT_EQ(jac.dintegrate_dq.cols(), 0);
    EXPECT_EQ(jac.dintegrate_dv.rows(), 0);
    EXPECT_EQ(jac.dintegrate_dv.cols(), 0);
  }
  EXPECT_EQ(jac.Phia().rows(), dimi);
  EXPECT_EQ(jac.Phia().cols(), dimv);
  EXPECT_EQ(jac.Phiq().rows(), dimi);
  EXPECT_EQ(jac.Phiq().cols(), dimv);
  EXPECT_EQ(jac.Phiv().rows(), dimi);
  EXPECT_EQ(jac.Phiv().cols(), dimv);
  EXPECT_EQ(jac.Phiu().rows(), dimi);
  EXPECT_EQ(jac.Phiu().cols(), dimu);
  jac.Phix().setRandom();
  EXPECT_TRUE(jac.Phix().leftCols(dimv).isApprox(jac.Phiq()));
  EXPECT_TRUE(jac.Phix().rightCols(dimv).isApprox(jac.Phiv()));
}


TEST_F(SplitStateConstraintJacobianTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.activateImpulse(0);
  test(robot, impulse_status);
}


TEST_F(SplitStateConstraintJacobianTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  auto impulse_status = robot.createImpulseStatus();
  test(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  test(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}