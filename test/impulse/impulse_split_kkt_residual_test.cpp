#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"


namespace idocp {

class SplitImpulseKKTResidualTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  static void testSize(const Robot& robot, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);

  virtual void TearDown() {
  }

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitImpulseKKTResidualTest::testSize(const Robot& robot, 
                                           const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = impulse_status.dimf();
  EXPECT_EQ(kkt_residual.dimf(), dimf);
  EXPECT_EQ(kkt_residual.splitKKTResidual().size(), 4*dimv+2*dimf);
  EXPECT_EQ(kkt_residual.Fq().size(), dimv);
  EXPECT_EQ(kkt_residual.Fv().size(), dimv);
  EXPECT_EQ(kkt_residual.V().size(), dimf);
  EXPECT_EQ(kkt_residual.lf().size(), dimf);
  EXPECT_EQ(kkt_residual.lq().size(), dimv);
  EXPECT_EQ(kkt_residual.lv().size(), dimv);
  EXPECT_EQ(kkt_residual.lx().size(), 2*dimv);
  EXPECT_EQ(kkt_residual.dimf(), dimf);
  EXPECT_EQ(kkt_residual.ldv.size(), dimv);
  kkt_residual.splitKKTResidual() = Eigen::VectorXd::Random(4*dimv+2*dimf);
  const Eigen::VectorXd ldv_ref = Eigen::VectorXd::Random(dimf);
  kkt_residual.ldv = ldv_ref;
  const Eigen::VectorXd Fq_ref = kkt_residual.splitKKTResidual().segment(0,             dimv);
  const Eigen::VectorXd Fv_ref = kkt_residual.splitKKTResidual().segment(dimv,          dimv);
  const Eigen::VectorXd Fx_ref = kkt_residual.splitKKTResidual().segment(0,           2*dimv);
  const Eigen::VectorXd V_ref  = kkt_residual.splitKKTResidual().segment(2*dimv,        dimf);
  const Eigen::VectorXd lf_ref = kkt_residual.splitKKTResidual().segment(2*dimv+dimf,   dimf);
  const Eigen::VectorXd lq_ref = kkt_residual.splitKKTResidual().segment(2*dimv+2*dimf, dimv);
  const Eigen::VectorXd lv_ref = kkt_residual.splitKKTResidual().segment(3*dimv+2*dimf, dimv);
  const Eigen::VectorXd lx_ref = kkt_residual.splitKKTResidual().segment(2*dimv+2*dimf, 2*dimv);
  EXPECT_TRUE(kkt_residual.Fq().isApprox(Fq_ref));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(Fv_ref));
  EXPECT_TRUE(kkt_residual.Fx().isApprox(Fx_ref));
  EXPECT_TRUE(kkt_residual.V().isApprox(V_ref));
  EXPECT_TRUE(kkt_residual.lf().isApprox(lf_ref));
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(kkt_residual.lx().isApprox(lx_ref));
  EXPECT_TRUE(kkt_residual.ldv.isApprox(ldv_ref));
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.splitKKTResidual().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitImpulseKKTResidualTest::testIsApprox(const Robot& robot, 
                                               const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.splitKKTResidual().setRandom();
  kkt_residual.ldv.setRandom();
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fq().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fq() = kkt_residual_ref.Fq();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fv().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fv() = kkt_residual_ref.Fv();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fx().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.Fx() = kkt_residual_ref.Fx();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lq().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lq() = kkt_residual_ref.lq();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lv().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lv() = kkt_residual_ref.lv();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lx().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.lx() = kkt_residual_ref.lx();
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.ldv.setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual.ldv = kkt_residual_ref.ldv;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  if (impulse_status.hasActiveImpulse()) {
    kkt_residual_ref.lf().setRandom();
    EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
    kkt_residual.lf() = kkt_residual_ref.lf();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
    kkt_residual_ref.V().setRandom();
    EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
    kkt_residual.V() = kkt_residual_ref.V();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  }
  else {
    kkt_residual_ref.lf().setRandom();
    kkt_residual_ref.V().setRandom();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  }
  kkt_residual_ref = kkt_residual;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.splitKKTResidual().setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitImpulseKKTResidualTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  auto impulse_status = robot.createImpulseStatus();
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.activateImpulse(0);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(SplitImpulseKKTResidualTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  auto impulse_status = robot.createImpulseStatus();
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}