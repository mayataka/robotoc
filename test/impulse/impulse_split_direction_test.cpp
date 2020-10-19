#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"


namespace idocp {

class ImpulseSplitDirectionTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(ImpulseSplitDirectionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true};
  ContactStatus contact_status(robot.max_point_contacts());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitDirection d(robot);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimc(), 0);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv());
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+2*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  d.setContactStatus(contact_status);
  EXPECT_EQ(d.dimf(), 3);
  EXPECT_EQ(d.dimc(), 3);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv()+2*3);
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+2*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), d.dimc());
  EXPECT_EQ(d.df().size(), d.dimf());
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimf = contact_status.dimf();
  const int dimc = dimf;
  const Eigen::VectorXd dlmd = split_direction.segment(                     0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(          robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(         2*robot.dimv(),          dimc);
  const Eigen::VectorXd dmu_position = split_direction.segment(2*robot.dimv(),          dimf);
  const Eigen::VectorXd dmu_velocity = split_direction.segment(2*robot.dimv()+dimf,     dimf);
  const Eigen::VectorXd df = split_direction.segment(     2*robot.dimv()+dimc,          dimf);
  const Eigen::VectorXd dq = split_direction.segment(2*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(3*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(2*robot.dimv()+dimc+dimf, 2*robot.dimv());
  const Eigen::VectorXd ddv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(robot.dimv());
  d.ddv = ddv;
  d.dbeta = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(ddv.isApprox(d.ddv));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  EXPECT_TRUE(d.dlmd().isZero());
  EXPECT_TRUE(d.dgmm().isZero());
  EXPECT_TRUE(d.dmu().isZero());
  EXPECT_TRUE(d.df().isZero());
  EXPECT_TRUE(d.dq().isZero());
  EXPECT_TRUE(d.dv().isZero());
  EXPECT_TRUE(d.dx().isZero());
  EXPECT_TRUE(d.ddv.isZero());
  EXPECT_TRUE(d.dbeta.isZero());
  const ImpulseSplitDirection d_random = ImpulseSplitDirection::Random(robot, contact_status);
  EXPECT_EQ(d_random.dimf(), 3);
  EXPECT_EQ(d_random.dimc(), 3);
  EXPECT_EQ(d_random.dimKKT(), 4*robot.dimv()+2*3);
  EXPECT_EQ(d_random.max_dimKKT(), 4*robot.dimv()+2*3);
  EXPECT_EQ(d_random.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), d.dimc());
  EXPECT_EQ(d_random.df().size(), d.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d_random.ddv.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
}


TEST_F(ImpulseSplitDirectionTest, floating_base_full_contacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active = {true, true, true, true};
  ContactStatus contact_status(robot.max_point_contacts());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitDirection d(robot);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimc(), 0);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv());
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  d.setContactStatus(contact_status);
  EXPECT_EQ(d.dimf(), 4*3);
  EXPECT_EQ(d.dimc(), 4*3);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), d.dimc());
  EXPECT_EQ(d.df().size(), d.dimf());
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimf = contact_status.dimf();
  const int dimc = dimf;
  const Eigen::VectorXd dlmd = split_direction.segment(                     0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(          robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(         2*robot.dimv(),          dimc);
  const Eigen::VectorXd dmu_position = split_direction.segment(2*robot.dimv(),          dimf);
  const Eigen::VectorXd dmu_velocity = split_direction.segment(2*robot.dimv()+dimf,     dimf);
  const Eigen::VectorXd df = split_direction.segment(     2*robot.dimv()+dimc,          dimf);
  const Eigen::VectorXd dq = split_direction.segment(2*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(3*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(2*robot.dimv()+dimc+dimf, 2*robot.dimv());
  const Eigen::VectorXd ddv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(robot.dimv());
  d.ddv = ddv;
  d.dbeta = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(ddv.isApprox(d.ddv));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  EXPECT_TRUE(d.dlmd().isZero());
  EXPECT_TRUE(d.dgmm().isZero());
  EXPECT_TRUE(d.dmu().isZero());
  EXPECT_TRUE(d.df().isZero());
  EXPECT_TRUE(d.dq().isZero());
  EXPECT_TRUE(d.dv().isZero());
  EXPECT_TRUE(d.dx().isZero());
  EXPECT_TRUE(d.ddv.isZero());
  EXPECT_TRUE(d.dbeta.isZero());
  const ImpulseSplitDirection d_random = ImpulseSplitDirection::Random(robot, contact_status);
  EXPECT_EQ(d.dimf(), 4*3);
  EXPECT_EQ(d.dimc(), 4*3);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d_random.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), d.dimc());
  EXPECT_EQ(d_random.df().size(), d.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d_random.ddv.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
}


TEST_F(ImpulseSplitDirectionTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  ContactStatus contact_status(robot.max_point_contacts());
  contact_status.setContactStatus(is_contact_active);
  ImpulseSplitDirection d(robot);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimc(), 0);
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv());
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), 0);
  EXPECT_EQ(d.df().size(), 0);
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  d.setContactStatus(contact_status);
  EXPECT_EQ(d.dimf(), contact_status.dimf());
  EXPECT_EQ(d.dimc(), contact_status.dimf());
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv()+2*contact_status.dimf());
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d.dlmd().size(), robot.dimv());
  EXPECT_EQ(d.dgmm().size(), robot.dimv());
  EXPECT_EQ(d.dmu().size(), d.dimc());
  EXPECT_EQ(d.df().size(), d.dimf());
  EXPECT_EQ(d.dq().size(), robot.dimv());
  EXPECT_EQ(d.dv().size(), robot.dimv());
  EXPECT_EQ(d.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d.ddv.size(), robot.dimv());
  EXPECT_EQ(d.dbeta.size(), robot.dimv());
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimf = contact_status.dimf();
  const int dimc = dimf;
  const Eigen::VectorXd dlmd = split_direction.segment(                     0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(          robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(         2*robot.dimv(),          dimc);
  const Eigen::VectorXd dmu_position = split_direction.segment(2*robot.dimv(),          dimf);
  const Eigen::VectorXd dmu_velocity = split_direction.segment(2*robot.dimv()+dimf,     dimf);
  const Eigen::VectorXd df = split_direction.segment(     2*robot.dimv()+dimc,          dimf);
  const Eigen::VectorXd dq = split_direction.segment(2*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(3*robot.dimv()+dimc+dimf,  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(2*robot.dimv()+dimc+dimf, 2*robot.dimv());
  const Eigen::VectorXd ddv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dbeta = Eigen::VectorXd::Random(robot.dimv());
  d.ddv = ddv;
  d.dbeta = dbeta;
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  EXPECT_TRUE(ddv.isApprox(d.ddv));
  EXPECT_TRUE(dbeta.isApprox(d.dbeta));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  EXPECT_TRUE(d.dlmd().isZero());
  EXPECT_TRUE(d.dgmm().isZero());
  EXPECT_TRUE(d.dmu().isZero());
  EXPECT_TRUE(d.df().isZero());
  EXPECT_TRUE(d.dq().isZero());
  EXPECT_TRUE(d.dv().isZero());
  EXPECT_TRUE(d.dx().isZero());
  EXPECT_TRUE(d.ddv.isZero());
  EXPECT_TRUE(d.dbeta.isZero());
  const ImpulseSplitDirection d_random = ImpulseSplitDirection::Random(robot, contact_status);
  EXPECT_EQ(d.dimf(), contact_status.dimf());
  EXPECT_EQ(d.dimc(), contact_status.dimf());
  EXPECT_EQ(d.dimKKT(), 4*robot.dimv()+2*contact_status.dimf());
  EXPECT_EQ(d.max_dimKKT(), 4*robot.dimv()+4*3+4*3);
  EXPECT_EQ(d_random.split_direction().size(), d.dimKKT());
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), d.dimc());
  EXPECT_EQ(d_random.df().size(), d.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.dx().size(), 2*robot.dimv());
  EXPECT_EQ(d_random.ddv.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}