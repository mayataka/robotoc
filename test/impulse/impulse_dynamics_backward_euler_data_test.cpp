#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler_data.hpp"


namespace idocp {

class ImpulseDynamicsForwardEulerDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void testSize(const Robot& robot, const ImpulseStatus& impulse_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseDynamicsForwardEulerDataTest::testSize(const Robot& robot, const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimp();
  ImpulseDynamicsBackwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  EXPECT_EQ(data.dImDdq.rows(), dimv);
  EXPECT_EQ(data.dImDdq.cols(), dimv);
  EXPECT_EQ(data.dImDddv.rows(), dimv);
  EXPECT_EQ(data.dImDddv.cols(), dimv);
  EXPECT_EQ(data.Minv.rows(), dimv);
  EXPECT_EQ(data.Minv.cols(), dimv);
  EXPECT_EQ(data.Qdvq.rows(), dimv);
  EXPECT_EQ(data.Qdvq.cols(), dimv);
  EXPECT_EQ(data.ImD.size(), dimv);
  EXPECT_EQ(data.Minv_ImD.size(), dimv);
  EXPECT_EQ(data.ldv.size(), dimv);
  EXPECT_EQ(data.Qdvf().rows(), dimv);
  EXPECT_EQ(data.Qdvf().cols(), dimf);

  const Eigen::MatrixXd Qdvf_ref = Eigen::MatrixXd::Random(dimv, dimf);
  data.Qdvf() = Qdvf_ref;
  EXPECT_TRUE(data.Qdvf().isApprox(Qdvf_ref));
}


TEST_F(ImpulseDynamicsForwardEulerDataTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ImpulseStatus impulse_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false});
  testSize(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testSize(robot, impulse_status);
}


TEST_F(ImpulseDynamicsForwardEulerDataTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ImpulseStatus impulse_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  impulse_status.setImpulseStatus({false, false, false, false});
  testSize(robot, impulse_status);
  std::random_device rnd;
  std::vector<bool> is_impulse_active;
  for (const auto frame : contact_frames) {
    is_impulse_active.push_back(rnd()%2==0);
  }
  impulse_status.setImpulseStatus(is_impulse_active);
  testSize(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}