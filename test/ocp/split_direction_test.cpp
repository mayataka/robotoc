#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

class SplitDirectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf, floating_base_urdf;
};


TEST_F(SplitDirectionTest, fixed_base) {
  Robot robot(fixed_base_urdf);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = 0;
  SplitDirection d(robot);
  EXPECT_EQ(d.dlmdgmm().size(), dimx);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.du().size(), dimu);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.du_passive.size(), 6);
  EXPECT_EQ(d.dnu_passive.size(), 6);
  EXPECT_EQ(d.daf().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimKKT(), 4*dimv+dimu);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction = split_direction;
  const Eigen::VectorXd dlmd = split_direction.segment(0,  dimv);
  const Eigen::VectorXd dgmm = split_direction.segment(dimv,  dimv);
  const Eigen::VectorXd du = split_direction.segment(2*dimv,  dimu);
  const Eigen::VectorXd dq = split_direction.segment(2*dimv+dimu,  dimv);
  const Eigen::VectorXd dv = split_direction.segment(2*dimv+dimu+dimv,  dimv);
  const Eigen::VectorXd dx = split_direction.segment(2*dimv+dimu, 2*dimv);
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(du.isApprox(d.du()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  const Eigen::VectorXd daf = Eigen::VectorXd::Random(dimv+dimf);
  d.daf() = daf;
  EXPECT_TRUE(daf.isApprox(d.daf()));
  EXPECT_TRUE(daf.head(dimv).isApprox(d.da()));
  EXPECT_TRUE(daf.tail(dimf).isApprox(d.df()));
  const Eigen::VectorXd dbetamu = Eigen::VectorXd::Random(dimv+dimf);
  d.dbetamu() = dbetamu;
  EXPECT_TRUE(dbetamu.isApprox(d.dbetamu()));
  EXPECT_TRUE(dbetamu.head(dimv).isApprox(d.dbetamu()));
  EXPECT_TRUE(dbetamu.tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.split_direction.isZero());
  const SplitDirection d_random = SplitDirection::Random(robot);
  EXPECT_EQ(d_random.dlmd().size(), dimv);
  EXPECT_EQ(d_random.dgmm().size(), dimv);
  EXPECT_EQ(d_random.du().size(), dimu);
  EXPECT_EQ(d_random.dq().size(), dimv);
  EXPECT_EQ(d_random.dv().size(), dimv);
  EXPECT_EQ(d_random.dx().size(), dimx);
  EXPECT_EQ(d_random.da().size(), dimv);
  EXPECT_EQ(d_random.dbeta().size(), dimv);
  EXPECT_EQ(d_random.du_passive.size(), 6);
  EXPECT_EQ(d_random.dnu_passive.size(), 6);
  EXPECT_EQ(d_random.daf().size(), dimv);
  EXPECT_EQ(d_random.dbetamu().size(), dimv);
  EXPECT_EQ(d_random.dmu().size(), 0);
  EXPECT_EQ(d_random.df().size(), 0);
  EXPECT_EQ(d_random.dimf(), 0);
  EXPECT_EQ(d_random.dimKKT(), 4*dimv+dimu);
  EXPECT_FALSE(d_random.split_direction.isZero());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.du().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.dx().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  EXPECT_FALSE(d_random.dbeta().isZero());
  EXPECT_TRUE(d_random.da().isApprox(d_random.daf().head(dimv)));
  EXPECT_TRUE(d_random.df().isApprox(d_random.daf().tail(dimf)));
  EXPECT_TRUE(d_random.dbeta().isApprox(d_random.dbetamu().head(dimv)));
  EXPECT_TRUE(d_random.dmu().isApprox(d_random.dbetamu().tail(dimf)));
}


TEST_F(SplitDirectionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  SplitDirection d(robot);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.du().size(), dimu);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.du_passive.size(), 6);
  EXPECT_EQ(d.dnu_passive.size(), 6);
  EXPECT_EQ(d.daf().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimKKT(), 4*dimv+dimu);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction = split_direction;
  const Eigen::VectorXd dlmd = split_direction.segment(0,  dimv);
  const Eigen::VectorXd dgmm = split_direction.segment(dimv,  dimv);
  const Eigen::VectorXd du = split_direction.segment(2*dimv,  dimu);
  const Eigen::VectorXd dq = split_direction.segment(2*dimv+dimu,  dimv);
  const Eigen::VectorXd dv = split_direction.segment(2*dimv+dimu+dimv,  dimv);
  const Eigen::VectorXd dx = split_direction.segment(2*dimv+dimu, 2*dimv);
  const Eigen::VectorXd da = Eigen::VectorXd::Random(dimv);
  d.da() = da;
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(dimv);
  d.dbeta() = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(du.isApprox(d.du()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(da.isApprox(d.daf()));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta()));
  EXPECT_TRUE(dbeta.isApprox(d.dbetamu()));
  EXPECT_TRUE(d.daf().head(dimv).isApprox(d.da()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  std::vector<bool> is_contact_active = {true};
  ContactStatus contact_status = ContactStatus(robot.max_point_contacts());
  contact_status.setContactStatus(is_contact_active);
  d.setContactStatus(contact_status);
  const int dimf = contact_status.dimf();
  EXPECT_EQ(d.daf().size(), dimv+dimf);
  EXPECT_EQ(d.df().size(), dimf);
  EXPECT_EQ(d.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d.dmu().size(), dimf);
  EXPECT_EQ(d.dimf(), dimf);
  const Eigen::VectorXd daf = Eigen::VectorXd::Random(dimv+dimf);
  d.daf() = daf;
  EXPECT_TRUE(daf.isApprox(d.daf()));
  EXPECT_TRUE(daf.head(dimv).isApprox(d.da()));
  EXPECT_TRUE(daf.tail(dimf).isApprox(d.df()));
  const Eigen::VectorXd dbetamu = Eigen::VectorXd::Random(dimv+dimf);
  d.dbetamu() = dbetamu;
  EXPECT_TRUE(dbetamu.isApprox(d.dbetamu()));
  EXPECT_TRUE(dbetamu.head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(dbetamu.tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.split_direction.isZero());
  const SplitDirection d_random = SplitDirection::Random(robot, contact_status);
  EXPECT_EQ(d_random.dlmd().size(), dimv);
  EXPECT_EQ(d_random.dgmm().size(), dimv);
  EXPECT_EQ(d_random.du().size(), dimu);
  EXPECT_EQ(d_random.dq().size(), dimv);
  EXPECT_EQ(d_random.dv().size(), dimv);
  EXPECT_EQ(d_random.dx().size(), dimx);
  EXPECT_EQ(d_random.da().size(), dimv);
  EXPECT_EQ(d_random.dbeta().size(), dimv);
  EXPECT_EQ(d_random.du_passive.size(), 6);
  EXPECT_EQ(d_random.dnu_passive.size(), 6);
  EXPECT_EQ(d_random.daf().size(), dimv+dimf);
  EXPECT_EQ(d_random.df().size(), dimf);
  EXPECT_EQ(d_random.dbetamu().size(), dimv+dimf);
  EXPECT_EQ(d_random.dmu().size(), dimf);
  EXPECT_EQ(d_random.dimf(), dimf);
  EXPECT_EQ(d_random.dimKKT(), 4*dimv+dimu);
  EXPECT_FALSE(d_random.split_direction.isZero());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.du().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.dx().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  EXPECT_FALSE(d_random.dbeta().isZero());
  EXPECT_FALSE(d_random.dmu().isZero());
  EXPECT_FALSE(d_random.df().isZero());
  EXPECT_TRUE(d_random.daf().head(dimv).isApprox(d_random.da()));
  EXPECT_TRUE(d_random.daf().tail(dimf).isApprox(d_random.df()));
  EXPECT_TRUE(d_random.dbetamu().head(dimv).isApprox(d_random.dbeta()));
  EXPECT_TRUE(d_random.dbetamu().tail(dimf).isApprox(d_random.dmu()));
}


TEST_F(SplitDirectionTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  SplitDirection d(robot);
  EXPECT_EQ(d.dlmd().size(), dimv);
  EXPECT_EQ(d.dgmm().size(), dimv);
  EXPECT_EQ(d.du().size(), dimu);
  EXPECT_EQ(d.dq().size(), dimv);
  EXPECT_EQ(d.dv().size(), dimv);
  EXPECT_EQ(d.dx().size(), dimx);
  EXPECT_EQ(d.da().size(), dimv);
  EXPECT_EQ(d.dbeta().size(), dimv);
  EXPECT_EQ(d.du_passive.size(), 6);
  EXPECT_EQ(d.dnu_passive.size(), 6);
  EXPECT_EQ(d.daf().size(), dimv);
  EXPECT_EQ(d.dbetamu().size(), dimv);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimKKT(), 4*dimv+dimu);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction = split_direction;
  const Eigen::VectorXd dlmd = split_direction.segment(0,  dimv);
  const Eigen::VectorXd dgmm = split_direction.segment(dimv,  dimv);
  const Eigen::VectorXd du = split_direction.segment(2*dimv,  dimu);
  const Eigen::VectorXd dq = split_direction.segment(2*dimv+dimu,  dimv);
  const Eigen::VectorXd dv = split_direction.segment(2*dimv+dimu+dimv,  dimv);
  const Eigen::VectorXd dx = split_direction.segment(2*dimv+dimu, 2*dimv);
  const Eigen::VectorXd da = Eigen::VectorXd::Random(dimv);
  d.da() = da;
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(dimv);
  d.dbeta() = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(du.isApprox(d.du()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(da.isApprox(d.daf()));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta()));
  EXPECT_TRUE(dbeta.isApprox(d.dbetamu()));
  EXPECT_TRUE(d.daf().head(dimv).isApprox(d.da()));
  EXPECT_TRUE(d.dbetamu().head(dimv).isApprox(d.dbeta()));
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  ContactStatus contact_status = ContactStatus(robot.max_point_contacts());
  contact_status.setContactStatus(is_contact_active);
  d.setContactStatus(contact_status);
  const int dimf = contact_status.dimf();
  EXPECT_EQ(d.dmu().size(), dimf);
  EXPECT_EQ(d.df().size(), dimf);
  EXPECT_EQ(d.dimf(), dimf);
  const Eigen::VectorXd daf = Eigen::VectorXd::Random(dimv+dimf);
  d.daf() = daf;
  EXPECT_TRUE(daf.isApprox(d.daf()));
  EXPECT_TRUE(daf.head(dimv).isApprox(d.da()));
  EXPECT_TRUE(daf.tail(dimf).isApprox(d.df()));
  const Eigen::VectorXd dbetamu = Eigen::VectorXd::Random(dimv+dimf);
  d.dbetamu() = dbetamu;
  EXPECT_TRUE(dbetamu.isApprox(d.dbetamu()));
  EXPECT_TRUE(dbetamu.head(dimv).isApprox(d.dbeta()));
  EXPECT_TRUE(dbetamu.tail(dimf).isApprox(d.dmu()));
  d.setZero();
  EXPECT_TRUE(d.split_direction.isZero());
  const SplitDirection d_random = SplitDirection::Random(robot, contact_status);
  EXPECT_EQ(d_random.dlmd().size(), dimv);
  EXPECT_EQ(d_random.dgmm().size(), dimv);
  EXPECT_EQ(d_random.du().size(), dimu);
  EXPECT_EQ(d_random.dq().size(), dimv);
  EXPECT_EQ(d_random.dv().size(), dimv);
  EXPECT_EQ(d_random.dx().size(), dimx);
  EXPECT_EQ(d_random.da().size(), dimv);
  EXPECT_EQ(d_random.dbeta().size(), dimv);
  EXPECT_EQ(d_random.du_passive.size(), 6);
  EXPECT_EQ(d_random.dnu_passive.size(), 6);
  EXPECT_EQ(d_random.dmu().size(), dimf);
  EXPECT_EQ(d_random.df().size(), dimf);
  EXPECT_EQ(d_random.dimf(), dimf);
  EXPECT_EQ(d_random.dimKKT(), 4*dimv+dimu);
  EXPECT_FALSE(d_random.split_direction.isZero());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.du().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.dx().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  EXPECT_FALSE(d_random.dbeta().isZero());
  EXPECT_FALSE(d_random.du_passive.isZero());
  EXPECT_FALSE(d_random.dnu_passive.isZero());
  EXPECT_TRUE(d_random.daf().head(dimv).isApprox(d_random.da()));
  EXPECT_TRUE(d_random.dbetamu().head(dimv).isApprox(d_random.dbeta()));
  if (contact_status.hasActiveContacts()) {
    EXPECT_FALSE(d_random.dmu().isZero());
    EXPECT_FALSE(d_random.df().isZero());
    EXPECT_TRUE(d_random.daf().tail(dimf).isApprox(d_random.df()));
    EXPECT_TRUE(d_random.dbetamu().tail(dimf).isApprox(d_random.dmu()));
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}