#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler_data.hpp"


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
  ImpulseDynamicsForwardEulerData data(robot);
  data.setImpulseStatus(impulse_status);
  EXPECT_EQ(data.dImDddv.rows(), dimv);
  EXPECT_EQ(data.dImDddv.cols(), dimv);
  EXPECT_EQ(data.dCddv().rows(), dimf);
  EXPECT_EQ(data.dCddv().cols(), dimv);

  EXPECT_EQ(data.dImDCdqv().rows(), dimv+dimf);
  EXPECT_EQ(data.dImDCdqv().cols(), dimx);

  EXPECT_EQ(data.dImDCdq().rows(), dimv+dimf);
  EXPECT_EQ(data.dImDCdq().cols(), dimv);

  EXPECT_EQ(data.dCdq().rows(), dimf);
  EXPECT_EQ(data.dCdq().cols(), dimv);

  EXPECT_EQ(data.dCdv().rows(), dimf);
  EXPECT_EQ(data.dCdv().cols(), dimv);

  EXPECT_EQ(data.MJtJinv().rows(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv().cols(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv_dImDCdqv().rows(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv_dImDCdqv().cols(), dimx);
  EXPECT_EQ(data.Qdvfqv().rows(), dimv+dimf);
  EXPECT_EQ(data.Qdvfqv().cols(), dimx);
  EXPECT_EQ(data.ImDC().size(), dimv+dimf);
  EXPECT_EQ(data.ImD().size(), dimv);
  EXPECT_EQ(data.C().size(), dimf);
  EXPECT_EQ(data.MJtJinv_ImDC().size(), dimv+dimf);
  EXPECT_EQ(data.ldvf().size(), dimv+dimf);
  EXPECT_EQ(data.ldv().size(), dimv);
  EXPECT_EQ(data.lf().size(), dimf);
  const Eigen::MatrixXd dImDddv_ref = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd dCddv_ref = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd dImDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd MJtJinv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
  const Eigen::MatrixXd MJtJinv_dImDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd Qdvfqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::VectorXd ImDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd MJtJinv_ImDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd ldvf_ref = Eigen::VectorXd::Random(dimv+dimf);
  data.dImDddv = dImDddv_ref;
  data.dCddv() = dCddv_ref;
  data.dImDCdqv() = dImDCdqv_ref;
  data.MJtJinv() = MJtJinv_ref;
  data.MJtJinv_dImDCdqv() = MJtJinv_dImDCdqv_ref;
  data.Qdvfqv() = Qdvfqv_ref;
  data.ImDC() = ImDC_ref;
  data.MJtJinv_ImDC() = MJtJinv_ImDC_ref;
  data.ldvf() = ldvf_ref;
  EXPECT_TRUE(data.dImDddv.isApprox(dImDddv_ref));
  EXPECT_TRUE(data.dCddv().isApprox(dCddv_ref));
  EXPECT_TRUE(data.dImDCdqv().isApprox(dImDCdqv_ref));
  EXPECT_TRUE(data.dImDCdq().isApprox(dImDCdqv_ref.leftCols(dimv)));
  EXPECT_TRUE(data.dImDdq().isApprox(dImDCdqv_ref.block(0, 0, dimv, dimv)));
  EXPECT_TRUE(data.dCdq().isApprox(dImDCdqv_ref.block(dimv, 0, dimf, dimv)));
  EXPECT_TRUE(data.dCdv().isApprox(dImDCdqv_ref.block(dimv, dimv, dimf, dimv)));
  EXPECT_TRUE(data.MJtJinv().isApprox(MJtJinv_ref));
  EXPECT_TRUE(data.MJtJinv_dImDCdqv().isApprox(MJtJinv_dImDCdqv_ref));
  EXPECT_TRUE(data.Qdvfqv().isApprox(Qdvfqv_ref));
  EXPECT_TRUE(data.ImDC().isApprox(ImDC_ref));
  EXPECT_TRUE(data.ImD().isApprox(ImDC_ref.head(dimv)));
  EXPECT_TRUE(data.C().isApprox(ImDC_ref.tail(dimf)));
  EXPECT_TRUE(data.MJtJinv_ImDC().isApprox(MJtJinv_ImDC_ref));
  EXPECT_TRUE(data.ldvf().isApprox(ldvf_ref));
  EXPECT_TRUE(data.ldv().isApprox(ldvf_ref.head(dimv)));
  EXPECT_TRUE(data.lf().isApprox(ldvf_ref.tail(dimf)));
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