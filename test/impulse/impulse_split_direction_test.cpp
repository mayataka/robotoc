#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

class ImpulseSplitDirectionTest : public ::testing::Test {
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
  static void testIsApprox(const Robot& robot, const ImpulseStatus& impulse_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseSplitDirectionTest::testSize(const Robot& robot, 
                                         const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimf();
  ImpulseSplitDirection d(robot);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.ddv().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.ddvf().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  d.setImpulseStatus(impulse_status);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.ddv().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.ddvf().size(), dimv+dimf);
  EXPECT_EQ(d.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d.df().size(), dimf);
  EXPECT_EQ(d.dmu().size(), dimf);
  EXPECT_EQ(d.dimf(), dimf);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(4*dimv);
  d.split_direction = split_direction;
  const Eigen::VectorXd dlmd = split_direction.segment(0,  dimv);
  const Eigen::VectorXd dgmm = split_direction.segment(dimv,  dimv);
  const Eigen::VectorXd dq = split_direction.segment(2*dimv,  dimv);
  const Eigen::VectorXd dv = split_direction.segment(2*dimv+dimv,  dimv);
  const Eigen::VectorXd dx = split_direction.segment(2*dimv, 2*dimv);
  const Eigen::VectorXd ddv = Eigen::VectorXd::Random(dimv);
  d.ddv() = ddv;
  const Eigen::VectorXd df = Eigen::VectorXd::Random(dimf);
  d.df() = df;
  const Eigen::VectorXd dmu = Eigen::VectorXd::Random(dimf);
  d.dmu() = dmu;
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(dimv);
  d.dbeta() = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(ddv.isApprox(d.ddv()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(d.ddvf().head(dimv).isApprox(d.ddv()));
  EXPECT_TRUE(d.ddvf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.split_direction.isZero());
  EXPECT_TRUE(d.dlmd().isZero());
  EXPECT_TRUE(d.dgmm().isZero());
  EXPECT_TRUE(d.dq().isZero());
  EXPECT_TRUE(d.dv().isZero());
  EXPECT_TRUE(d.dx().isZero());
  EXPECT_TRUE(d.ddv().isZero());
  EXPECT_TRUE(d.df().isZero());
  EXPECT_TRUE(d.ddvf().isZero());
  EXPECT_TRUE(d.dbeta().isZero());
  EXPECT_TRUE(d.dmu().isZero());
  EXPECT_TRUE(d.dbetamu().isZero());
  d.setRandom();
  EXPECT_FALSE(d.split_direction.isZero());
  EXPECT_FALSE(d.dlmd().isZero());
  EXPECT_FALSE(d.dgmm().isZero());
  EXPECT_FALSE(d.dq().isZero());
  EXPECT_FALSE(d.dv().isZero());
  EXPECT_FALSE(d.dx().isZero());
  EXPECT_FALSE(d.ddv().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d.df().isZero());
  EXPECT_FALSE(d.ddvf().isZero());
  EXPECT_FALSE(d.dbeta().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d.dmu().isZero());
  EXPECT_FALSE(d.dbetamu().isZero());
  EXPECT_TRUE(d.ddvf().head(dimv).isApprox(d.ddv()));
  EXPECT_TRUE(d.ddvf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  const ImpulseSplitDirection d_random = ImpulseSplitDirection::Random(robot, impulse_status);
  EXPECT_EQ(d_random.dlmd().size(), dimv);
  EXPECT_EQ(d_random.dgmm().size(), dimv);
  EXPECT_EQ(d_random.dq().size(), dimv);
  EXPECT_EQ(d_random.dv().size(), dimv);
  EXPECT_EQ(d_random.dx().size(), dimx);
  EXPECT_EQ(d_random.ddv().size(), dimv);
  EXPECT_EQ(d_random.dbeta().size(), dimv);
  EXPECT_EQ(d_random.ddvf().size(), dimv+dimf);
  EXPECT_EQ(d_random.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d_random.df().size(), dimf);
  EXPECT_EQ(d_random.dmu().size(), dimf);
  EXPECT_EQ(d_random.dimf(), dimf);
  EXPECT_FALSE(d_random.split_direction.isZero());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.dx().isZero());
  EXPECT_FALSE(d_random.ddv().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d_random.df().isZero());
  EXPECT_FALSE(d_random.ddvf().isZero());
  EXPECT_FALSE(d_random.dbeta().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d_random.dmu().isZero());
  EXPECT_FALSE(d_random.dbetamu().isZero());
  EXPECT_TRUE(d_random.ddvf().head(dimv).isApprox(d_random.ddv()));
  EXPECT_TRUE(d_random.ddvf().tail(dimf).isApprox(d_random.df()));
  EXPECT_TRUE(d_random.dbetamu().head(dimv).isApprox(d_random.dbeta()));
  EXPECT_TRUE(d_random.dbetamu().tail(dimf).isApprox(d_random.dmu()));
}


void ImpulseSplitDirectionTest::testIsApprox(const Robot& robot, 
                                             const ImpulseStatus& impulse_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = impulse_status.dimf();
  ImpulseSplitDirection d(robot);
  d.setRandom(impulse_status);
  EXPECT_FALSE(d.split_direction.isZero());
  EXPECT_FALSE(d.dlmd().isZero());
  EXPECT_FALSE(d.dgmm().isZero());
  EXPECT_FALSE(d.dq().isZero());
  EXPECT_FALSE(d.dv().isZero());
  EXPECT_FALSE(d.dx().isZero());
  EXPECT_FALSE(d.ddv().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d.df().isZero());
  EXPECT_FALSE(d.ddvf().isZero());
  EXPECT_FALSE(d.dbeta().isZero());
  if (dimf > 0) 
  EXPECT_FALSE(d.dmu().isZero());
  EXPECT_FALSE(d.dbetamu().isZero());
  EXPECT_TRUE(d.ddvf().head(dimv).isApprox(d.ddv()));
  EXPECT_TRUE(d.ddvf().tail(dimf).isApprox(d.df()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(d.dbetamu().tail(dimf).isApprox(d.dmu()));
  ImpulseSplitDirection d_ref = d;
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dlmd().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dlmd() = d.dlmd();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dgmm().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dgmm() = d.dgmm();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dq().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dq() = d.dq();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dv().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dv() = d.dv();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.ddv().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.ddv() = d.ddv();
  EXPECT_TRUE(d.isApprox(d_ref));
  d_ref.dbeta().setRandom();
  EXPECT_FALSE(d.isApprox(d_ref));
  d_ref.dbeta() = d.dbeta();
  EXPECT_TRUE(d.isApprox(d_ref));
  if (impulse_status.hasActiveImpulse()) {
    d_ref.df().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.df() = d.df();
    EXPECT_TRUE(d.isApprox(d_ref));
    d_ref.dmu().setRandom();
    EXPECT_FALSE(d.isApprox(d_ref));
    d_ref.dmu() = d.dmu();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
  else {
    d_ref.df().setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
    d_ref.dmu().setRandom();
    EXPECT_TRUE(d.isApprox(d_ref));
  }
}


TEST_F(ImpulseSplitDirectionTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false};
  ImpulseStatus impulse_status = ImpulseStatus(robot.maxPointContacts());
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  impulse_status.activateImpulse(0);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
}


TEST_F(ImpulseSplitDirectionTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ImpulseStatus impulse_status = ImpulseStatus(robot.maxPointContacts());
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, impulse_status);
  testIsApprox(robot, impulse_status);
  is_contact_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
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