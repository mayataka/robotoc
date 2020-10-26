#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/kkt_residual.hpp"


namespace idocp {

class KKTResidualTest : public ::testing::Test {
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


TEST_F(KKTResidualTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active = {rnd()%2==0};
  contact_status.setContactStatus(is_contact_active);
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimu);
  EXPECT_EQ(residual.KKT_residual.size(), 4*dimv+dimu);
  EXPECT_EQ(residual.Fq().size(), dimv);
  EXPECT_EQ(residual.Fv().size(), dimv);
  EXPECT_EQ(residual.lu().size(), dimu);
  EXPECT_EQ(residual.lq().size(), dimv);
  EXPECT_EQ(residual.lv().size(), dimv);
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_EQ(residual.C().size(), contact_status.dimf());
  EXPECT_EQ(residual.lf().size(), contact_status.dimf());
  EXPECT_EQ(residual.dimf(), contact_status.dimf());
  EXPECT_EQ(residual.la.size(), dimv);
  EXPECT_EQ(residual.ID.size(), dimv);
  EXPECT_EQ(residual.lu_passive.size(), 0);
  EXPECT_EQ(residual.C_passive.size(), 0);
  residual.KKT_residual = Eigen::VectorXd::Random(residual.dimKKT());
  const Eigen::VectorXd Fq = residual.KKT_residual.segment(0, dimv);
  const Eigen::VectorXd Fv = residual.KKT_residual.segment(dimv, dimv);
  const Eigen::VectorXd lu = residual.KKT_residual.segment(2*dimv, dimu);
  const Eigen::VectorXd lq = residual.KKT_residual.segment(2*dimv+dimu, dimv);
  const Eigen::VectorXd lv = residual.KKT_residual.segment(3*dimv+dimu, dimv);
  const Eigen::VectorXd lx = residual.KKT_residual.segment(2*dimv+dimu, 2*dimv);
  EXPECT_TRUE(residual.Fq().isApprox(Fq));
  EXPECT_TRUE(residual.Fv().isApprox(Fv));
  EXPECT_TRUE(residual.lu().isApprox(lu));
  EXPECT_TRUE(residual.lq().isApprox(lq));
  EXPECT_TRUE(residual.lv().isApprox(lv));
  EXPECT_TRUE(residual.lx().isApprox(lx));
}


TEST_F(KKTResidualTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  EXPECT_EQ(residual.dimf(), dimf);
  EXPECT_EQ(residual.dimKKT(), 4*dimv+dimu);
  EXPECT_EQ(residual.KKT_residual.size(), 4*dimv+dimu);
  EXPECT_EQ(residual.Fq().size(), dimv);
  EXPECT_EQ(residual.Fv().size(), dimv);
  EXPECT_EQ(residual.lu().size(), dimu);
  EXPECT_EQ(residual.lq().size(), dimv);
  EXPECT_EQ(residual.lv().size(), dimv);
  EXPECT_EQ(residual.lx().size(), 2*dimv);
  EXPECT_EQ(residual.C().size(), contact_status.dimf());
  EXPECT_EQ(residual.lf().size(), contact_status.dimf());
  EXPECT_EQ(residual.dimf(), contact_status.dimf());
  EXPECT_EQ(residual.la.size(), dimv);
  EXPECT_EQ(residual.ID.size(), dimv);
  EXPECT_EQ(residual.lu_passive.size(), 6);
  EXPECT_EQ(residual.C_passive.size(), 6);
  residual.KKT_residual = Eigen::VectorXd::Random(residual.dimKKT());
  const Eigen::VectorXd Fq = residual.KKT_residual.segment(0, dimv);
  const Eigen::VectorXd Fv = residual.KKT_residual.segment(dimv, dimv);
  const Eigen::VectorXd lu = residual.KKT_residual.segment(2*dimv, dimu);
  const Eigen::VectorXd lq = residual.KKT_residual.segment(2*dimv+dimu, dimv);
  const Eigen::VectorXd lv = residual.KKT_residual.segment(3*dimv+dimu, dimv);
  const Eigen::VectorXd lx = residual.KKT_residual.segment(2*dimv+dimu, 2*dimv);
  EXPECT_TRUE(residual.Fq().isApprox(Fq));
  EXPECT_TRUE(residual.Fv().isApprox(Fv));
  EXPECT_TRUE(residual.lu().isApprox(lu));
  EXPECT_TRUE(residual.lq().isApprox(lq));
  EXPECT_TRUE(residual.lv().isApprox(lv));
  EXPECT_TRUE(residual.lx().isApprox(lx));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}