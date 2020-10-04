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
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
};


TEST_F(SplitDirectionTest, fixed_base) {
  Robot robot(fixed_base_urdf_);
  std::random_device rnd;
  SplitDirection d(robot);
  EXPECT_EQ(d.dimf(), 0);
  EXPECT_EQ(d.dimc(), 0);
  EXPECT_EQ(d.dimKKT(), 5*robot.dimv());
  EXPECT_EQ(d.max_dimKKT(), 5*robot.dimv());
  d.setContactStatus(robot);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimc = robot.dim_passive() + robot.dimf();
  const Eigen::VectorXd dlmd = split_direction.segment(                             0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(                  robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(                 2*robot.dimv(),          dimc);
  const Eigen::VectorXd da = split_direction.segment(             2*robot.dimv()+dimc,  robot.dimv());
  const Eigen::VectorXd df = split_direction.segment(             3*robot.dimv()+dimc,  robot.dimf());
  const Eigen::VectorXd dq = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(4*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(), 2*robot.dimv());
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  const SplitDirection d_random = SplitDirection::Random(robot);
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(d_random.da().size(), robot.dimv());
  EXPECT_EQ(d_random.df().size(), robot.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.du.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.du.isZero());
  EXPECT_FALSE(d_random.dbeta.isZero());
  EXPECT_EQ(d_random.dimf(), 0);
  EXPECT_EQ(d_random.dimc(), 0);
  EXPECT_EQ(d_random.dimKKT(), 5*robot.dimv());
  EXPECT_EQ(d_random.max_dimKKT(), 5*robot.dimv());
}


TEST_F(SplitDirectionTest, fixed_base_contact) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status = {true};
  robot.setContactStatus(contact_status);
  SplitDirection d(robot);
  EXPECT_EQ(d.dimf(), robot.dimf());
  EXPECT_EQ(d.dimc(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(d.dimKKT(), 5*robot.dimv()+2*robot.dimf());
  EXPECT_EQ(d.max_dimKKT(), 5*robot.dimv()+2*robot.max_dimf());
  d.setContactStatus(robot);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimc = robot.dim_passive() + robot.dimf();
  const Eigen::VectorXd dlmd = split_direction.segment(                             0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(                  robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(                 2*robot.dimv(),          dimc);
  const Eigen::VectorXd da = split_direction.segment(             2*robot.dimv()+dimc,  robot.dimv());
  const Eigen::VectorXd df = split_direction.segment(             3*robot.dimv()+dimc,  robot.dimf());
  const Eigen::VectorXd dq = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(4*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(), 2*robot.dimv());
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  const SplitDirection d_random = SplitDirection::Random(robot);
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(d_random.da().size(), robot.dimv());
  EXPECT_EQ(d_random.df().size(), robot.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.du.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.dmu().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  EXPECT_FALSE(d_random.df().isZero());
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.du.isZero());
  EXPECT_FALSE(d_random.dbeta.isZero());
  EXPECT_EQ(d_random.dimf(), 3);
  EXPECT_EQ(d_random.dimc(), 3);
  EXPECT_EQ(d_random.dimKKT(), 5*robot.dimv()+2*robot.dimf());
  EXPECT_EQ(d_random.max_dimKKT(), 5*robot.dimv()+2*robot.max_dimf());
}


TEST_F(SplitDirectionTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> contact_status;
  for (const auto frame : contact_frames) {
    contact_status.push_back(rnd()%2==0);
  }
  robot.setContactStatus(contact_status);
  SplitDirection d(robot);
  EXPECT_EQ(d.dimKKT(), 5*robot.dimv()+robot.dim_passive()+2*robot.dimf());
  EXPECT_EQ(d.max_dimKKT(), 5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf());
  d.setContactStatus(robot);
  const Eigen::VectorXd split_direction = Eigen::VectorXd::Random(d.dimKKT());
  d.split_direction() = split_direction;
  const int dimc = robot.dim_passive() + robot.dimf();
  const Eigen::VectorXd dlmd = split_direction.segment(                             0,  robot.dimv());
  const Eigen::VectorXd dgmm = split_direction.segment(                  robot.dimv(),  robot.dimv());
  const Eigen::VectorXd dmu = split_direction.segment(                 2*robot.dimv(),          dimc);
  const Eigen::VectorXd da = split_direction.segment(             2*robot.dimv()+dimc,  robot.dimv());
  const Eigen::VectorXd df = split_direction.segment(             3*robot.dimv()+dimc,  robot.dimf());
  const Eigen::VectorXd dq = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dv = split_direction.segment(4*robot.dimv()+dimc+robot.dimf(),  robot.dimv());
  const Eigen::VectorXd dx = split_direction.segment(3*robot.dimv()+dimc+robot.dimf(), 2*robot.dimv());
  EXPECT_TRUE(dlmd.isApprox(d.dlmd()));
  EXPECT_TRUE(dgmm.isApprox(d.dgmm()));
  EXPECT_TRUE(dmu.isApprox(d.dmu()));
  EXPECT_TRUE(da.isApprox(d.da()));
  EXPECT_TRUE(df.isApprox(d.df()));
  EXPECT_TRUE(dq.isApprox(d.dq()));
  EXPECT_TRUE(dv.isApprox(d.dv()));
  EXPECT_TRUE(dx.isApprox(d.dx()));
  d.setZero();
  EXPECT_TRUE(d.split_direction().isZero());
  const SplitDirection d_random = SplitDirection::Random(robot);
  EXPECT_EQ(d_random.dlmd().size(), robot.dimv());
  EXPECT_EQ(d_random.dgmm().size(), robot.dimv());
  EXPECT_EQ(d_random.dmu().size(), robot.dim_passive()+robot.dimf());
  EXPECT_EQ(d_random.da().size(), robot.dimv());
  EXPECT_EQ(d_random.df().size(), robot.dimf());
  EXPECT_EQ(d_random.dq().size(), robot.dimv());
  EXPECT_EQ(d_random.dv().size(), robot.dimv());
  EXPECT_EQ(d_random.du.size(), robot.dimv());
  EXPECT_EQ(d_random.dbeta.size(), robot.dimv());
  EXPECT_FALSE(d_random.dlmd().isZero());
  EXPECT_FALSE(d_random.dgmm().isZero());
  EXPECT_FALSE(d_random.dmu().isZero());
  EXPECT_FALSE(d_random.da().isZero());
  if (robot.dimf() > 0) {
    EXPECT_FALSE(d_random.df().isZero());
  }
  EXPECT_FALSE(d_random.dq().isZero());
  EXPECT_FALSE(d_random.dv().isZero());
  EXPECT_FALSE(d_random.du.isZero());
  EXPECT_FALSE(d_random.dbeta.isZero());
  EXPECT_EQ(d_random.dimf(), robot.dimf());
  EXPECT_EQ(d_random.dimc(), 6+robot.dimf());
  EXPECT_EQ(d_random.dimKKT(), 5*robot.dimv()+6+2*robot.dimf());
  EXPECT_EQ(d_random.max_dimKKT(), 5*robot.dimv()+6+2*robot.max_dimf());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}