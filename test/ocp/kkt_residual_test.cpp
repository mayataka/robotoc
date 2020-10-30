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


TEST_F(KKTResidualTest, fixedBase) {
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
  EXPECT_EQ(residual.lf().size(), contact_status.dimf());
  EXPECT_EQ(residual.dimf(), contact_status.dimf());
  EXPECT_EQ(residual.la.size(), dimv);
  EXPECT_EQ(residual.lu_passive.size(), 6);
  residual.KKT_residual = Eigen::VectorXd::Random(residual.dimKKT());
  const Eigen::VectorXd la_ref = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lf_ref = Eigen::VectorXd::Random(contact_status.dimf());
  residual.la = la_ref;
  residual.lf() = lf_ref;
  const Eigen::VectorXd Fq_ref = residual.KKT_residual.segment(0, dimv);
  const Eigen::VectorXd Fv_ref = residual.KKT_residual.segment(dimv, dimv);
  const Eigen::VectorXd lu_ref = residual.KKT_residual.segment(2*dimv, dimu);
  const Eigen::VectorXd lq_ref = residual.KKT_residual.segment(2*dimv+dimu, dimv);
  const Eigen::VectorXd lv_ref = residual.KKT_residual.segment(3*dimv+dimu, dimv);
  const Eigen::VectorXd lx_ref = residual.KKT_residual.segment(2*dimv+dimu, 2*dimv);
  EXPECT_TRUE(residual.Fq().isApprox(Fq_ref));
  EXPECT_TRUE(residual.Fv().isApprox(Fv_ref));
  EXPECT_TRUE(residual.lu().isApprox(lu_ref));
  EXPECT_TRUE(residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(residual.lx().isApprox(lx_ref));
  EXPECT_TRUE(residual.la.isApprox(la_ref));
  EXPECT_TRUE(residual.lf().isApprox(lf_ref));
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.KKT_residual.setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}


TEST_F(KKTResidualTest, isApproxFixedBaseWithoutContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus({false});
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  residual.KKT_residual.setRandom();
  residual.la.setRandom();
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lf().setRandom();
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lu_passive.setRandom();
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.KKT_residual.setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}


TEST_F(KKTResidualTest, isApproxFixedBaseWithContacts) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus({true});
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  residual.KKT_residual.setRandom();
  residual.la.setRandom();
  residual.lf().setRandom();
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lu_passive.setRandom();
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lf().setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}


TEST_F(KKTResidualTest, floatingBase) {
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
  EXPECT_EQ(residual.lf().size(), contact_status.dimf());
  EXPECT_EQ(residual.dimf(), contact_status.dimf());
  EXPECT_EQ(residual.la.size(), dimv);
  EXPECT_EQ(residual.lu_passive.size(), 6);
  residual.KKT_residual = Eigen::VectorXd::Random(residual.dimKKT());
  const Eigen::VectorXd la_ref = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd lf_ref = Eigen::VectorXd::Random(contact_status.dimf());
  residual.la = la_ref;
  residual.lf() = lf_ref;
  const Eigen::VectorXd Fq_ref = residual.KKT_residual.segment(0, dimv);
  const Eigen::VectorXd Fv_ref = residual.KKT_residual.segment(dimv, dimv);
  const Eigen::VectorXd lu_ref = residual.KKT_residual.segment(2*dimv, dimu);
  const Eigen::VectorXd lq_ref = residual.KKT_residual.segment(2*dimv+dimu, dimv);
  const Eigen::VectorXd lv_ref = residual.KKT_residual.segment(3*dimv+dimu, dimv);
  const Eigen::VectorXd lx_ref = residual.KKT_residual.segment(2*dimv+dimu, 2*dimv);
  EXPECT_TRUE(residual.Fq().isApprox(Fq_ref));
  EXPECT_TRUE(residual.Fv().isApprox(Fv_ref));
  EXPECT_TRUE(residual.lu().isApprox(lu_ref));
  EXPECT_TRUE(residual.lq().isApprox(lq_ref));
  EXPECT_TRUE(residual.lv().isApprox(lv_ref));
  EXPECT_TRUE(residual.lx().isApprox(lx_ref));
  EXPECT_TRUE(residual.la.isApprox(la_ref));
  EXPECT_TRUE(residual.lf().isApprox(lf_ref));
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.KKT_residual.setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}


TEST_F(KKTResidualTest, isApproxFloatingBaseWithoutContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus({false, false, false, false});
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  residual.KKT_residual.setRandom();
  residual.la.setRandom();
  residual.lu_passive.setRandom();
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lf().setRandom();
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lu_passive.setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}


TEST_F(KKTResidualTest, isApproxFloatingBaseWithContacts) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus(is_contact_active);
  KKTResidual residual(robot);
  residual.setContactStatus(contact_status);
  residual.KKT_residual.setRandom();
  residual.la.setRandom();
  residual.lf().setRandom();
  residual.lu_passive.setRandom();
  KKTResidual kkt_residual_ref = residual;
  EXPECT_TRUE(residual.isApprox(kkt_residual_ref));
  kkt_residual_ref.lf().setRandom();
  EXPECT_FALSE(residual.isApprox(kkt_residual_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}