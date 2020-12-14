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
  const int dimf = impulse_status.dimp();
  EXPECT_EQ(kkt_residual.dimf(), dimf);
  EXPECT_EQ(kkt_residual.dimKKT(), 4*dimv);
  EXPECT_EQ(kkt_residual.KKT_residual.size(), 4*dimv);
  EXPECT_EQ(kkt_residual.Fq().size(), dimv);
  EXPECT_EQ(kkt_residual.Fv().size(), dimv);
  EXPECT_EQ(kkt_residual.P().size(), dimf);
  EXPECT_EQ(kkt_residual.V().size(), dimf);
  EXPECT_EQ(kkt_residual.lq().size(), dimv);
  EXPECT_EQ(kkt_residual.lv().size(), dimv);
  EXPECT_EQ(kkt_residual.lx().size(), 2*dimv);
  EXPECT_EQ(kkt_residual.lf().size(), dimf);
  EXPECT_EQ(kkt_residual.dimf(), dimf);
  EXPECT_EQ(kkt_residual.ldv.size(), dimv);
  kkt_residual.KKT_residual = Eigen::VectorXd::Random(kkt_residual.dimKKT());
  const Eigen::VectorXd P_ref = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd V_ref = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd ldv_ref = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lf_ref = Eigen::VectorXd::Random(dimf);
  kkt_residual.P() = P_ref;
  kkt_residual.V() = V_ref;
  kkt_residual.ldv = ldv_ref;
  kkt_residual.lf() = lf_ref;
  const Eigen::VectorXd Fq_ref = kkt_residual.KKT_residual.segment(0, dimv);
  const Eigen::VectorXd Fv_ref = kkt_residual.KKT_residual.segment(dimv, dimv);
  const Eigen::VectorXd lq_ref = kkt_residual.KKT_residual.segment(2*dimv, dimv);
  const Eigen::VectorXd lv_ref = kkt_residual.KKT_residual.segment(3*dimv, dimv);
  const Eigen::VectorXd lx_ref = kkt_residual.KKT_residual.segment(2*dimv, 2*dimv);
  EXPECT_TRUE(kkt_residual.Fq().isApprox(Fq_ref));
  EXPECT_TRUE(kkt_residual.Fv().isApprox(Fv_ref));
  EXPECT_TRUE(kkt_residual.P().isApprox(P_ref));
  EXPECT_TRUE(kkt_residual.V().isApprox(V_ref));
  EXPECT_TRUE(kkt_residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(kkt_residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(kkt_residual.lx().isApprox(lx_ref));
  EXPECT_TRUE(kkt_residual.ldv.isApprox(ldv_ref));
  EXPECT_TRUE(kkt_residual.lf().isApprox(lf_ref));
  ImpulseSplitKKTResidual kkt_residual_ref = kkt_residual;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.KKT_residual.setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
}


void SplitImpulseKKTResidualTest::testIsApprox(const Robot& robot, 
                                               const ImpulseStatus& impulse_status) {
  ImpulseSplitKKTResidual kkt_residual(robot);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.KKT_residual.setRandom();
  kkt_residual.P().setRandom();
  kkt_residual.V().setRandom();
  kkt_residual.ldv.setRandom();
  kkt_residual.lf().setRandom();
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

    kkt_residual_ref.P().setRandom();
    EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
    kkt_residual.P() = kkt_residual_ref.P();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));

    kkt_residual_ref.V().setRandom();
    EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
    kkt_residual.V() = kkt_residual_ref.V();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  }
  else {
    kkt_residual_ref.lf().setRandom();
    kkt_residual_ref.P().setRandom();
    kkt_residual_ref.V().setRandom();
    EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  }
  kkt_residual_ref = kkt_residual;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.KKT_residual.setRandom();
  EXPECT_FALSE(kkt_residual.isApprox(kkt_residual_ref));
}


TEST_F(SplitImpulseKKTResidualTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  ImpulseStatus impulse_status(contact_frames.size());
  impulse_status.setImpulseStatus({false});
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(SplitImpulseKKTResidualTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ImpulseStatus impulse_status(contact_frames.size());
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  is_contact_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}