#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/contact_dynamics_data.hpp"


namespace idocp {

class ContactDynamicsDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void testSize(const Robot& robot, const ContactStatus& contact_status);

  std::string fixed_base_urdf, floating_base_urdf;
};


void ContactDynamicsDataTest::testSize(const Robot& robot, const ContactStatus& contact_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dim_passive = robot.dim_passive();
  ContactDynamicsData data(robot);
  data.setContactStatus(contact_status);
  EXPECT_EQ(data.dIDda.rows(), dimv);
  EXPECT_EQ(data.dIDda.cols(), dimv);
  EXPECT_EQ(data.u_passive.size(), 6);
  EXPECT_EQ(data.dCda().rows(), dimf);
  EXPECT_EQ(data.dCda().cols(), dimv);
  EXPECT_EQ(data.dIDCdqv().rows(), dimv+dimf);
  EXPECT_EQ(data.dIDCdqv().cols(), dimx);
  EXPECT_EQ(data.dIDdq().rows(), dimv);
  EXPECT_EQ(data.dIDdq().cols(), dimv);
  EXPECT_EQ(data.dIDdv().rows(), dimv);
  EXPECT_EQ(data.dIDdv().cols(), dimv);
  EXPECT_EQ(data.dCdq().rows(), dimf);
  EXPECT_EQ(data.dCdq().cols(), dimv);
  EXPECT_EQ(data.dCdv().rows(), dimf);
  EXPECT_EQ(data.dCdv().cols(), dimv);
  EXPECT_EQ(data.MJtJinv().rows(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv().cols(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv_dIDCdqv().rows(), dimv+dimf);
  EXPECT_EQ(data.MJtJinv_dIDCdqv().cols(), dimx);
  EXPECT_EQ(data.Qafqv().rows(), dimv+dimf);
  EXPECT_EQ(data.Qafqv().cols(), dimx);
  EXPECT_EQ(data.Qafu_full().rows(), dimv+dimf);
  EXPECT_EQ(data.Qafu_full().cols(), dimv);
  EXPECT_EQ(data.Qafu_passive().rows(), dimv+dimf);
  EXPECT_EQ(data.Qafu_passive().cols(), dim_passive);
  EXPECT_EQ(data.Qafu().rows(), dimv+dimf);
  EXPECT_EQ(data.Qafu().cols(), dimu);
  EXPECT_EQ(data.IDC().size(), dimv+dimf);
  EXPECT_EQ(data.ID_full().size(), dimv);
  EXPECT_EQ(data.ID_passive().size(), dim_passive);
  EXPECT_EQ(data.ID().size(), dimu);
  EXPECT_EQ(data.C().size(), dimf);
  EXPECT_EQ(data.MJtJinv_IDC().size(), dimv+dimf);
  EXPECT_EQ(data.laf().size(), dimv+dimf);
  EXPECT_EQ(data.la().size(), dimv);
  EXPECT_EQ(data.lf().size(), dimf);
  const Eigen::MatrixXd dIDda_ref = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::VectorXd u_passive_ref = Eigen::VectorXd::Random(6);
  const Eigen::MatrixXd dCda_ref = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd dIDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd MJtJinv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
  const Eigen::MatrixXd MJtJinv_dIDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd Qafqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd Qafu_full_ref = Eigen::MatrixXd::Random(dimv+dimf, dimv);
  const Eigen::VectorXd IDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd MJtJinv_IDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd laf_ref = Eigen::VectorXd::Random(dimv+dimf);
  data.dIDda = dIDda_ref;
  data.u_passive = u_passive_ref;
  data.dCda() = dCda_ref;
  data.dIDCdqv() = dIDCdqv_ref;
  data.MJtJinv() = MJtJinv_ref;
  data.MJtJinv_dIDCdqv() = MJtJinv_dIDCdqv_ref;
  data.Qafqv() = Qafqv_ref;
  data.Qafu_full() = Qafu_full_ref;
  data.IDC() = IDC_ref;
  data.MJtJinv_IDC() = MJtJinv_IDC_ref;
  data.laf() = laf_ref;
  EXPECT_TRUE(data.dIDda.isApprox(dIDda_ref));
  EXPECT_TRUE(data.u_passive.isApprox(u_passive_ref));
  EXPECT_TRUE(data.dCda().isApprox(dCda_ref));
  EXPECT_TRUE(data.dIDCdqv().isApprox(dIDCdqv_ref));
  EXPECT_TRUE(data.dIDdq().isApprox(dIDCdqv_ref.topLeftCorner(dimv, dimv)));
  EXPECT_TRUE(data.dIDdv().isApprox(dIDCdqv_ref.topRightCorner(dimv, dimv)));
  EXPECT_TRUE(data.dCdq().isApprox(dIDCdqv_ref.bottomLeftCorner(dimf, dimv)));
  EXPECT_TRUE(data.dCdv().isApprox(dIDCdqv_ref.bottomRightCorner(dimf, dimv)));
  EXPECT_TRUE(data.MJtJinv().isApprox(MJtJinv_ref));
  EXPECT_TRUE(data.MJtJinv_dIDCdqv().isApprox(MJtJinv_dIDCdqv_ref));
  EXPECT_TRUE(data.Qafqv().isApprox(Qafqv_ref));
  EXPECT_TRUE(data.Qafu_full().isApprox(Qafu_full_ref));
  EXPECT_TRUE(data.Qafu_passive().isApprox(Qafu_full_ref.leftCols(dim_passive)));
  EXPECT_TRUE(data.Qafu().isApprox(Qafu_full_ref.rightCols(dimu)));
  EXPECT_TRUE(data.IDC().isApprox(IDC_ref));
  EXPECT_TRUE(data.ID_full().isApprox(IDC_ref.head(dimv)));
  EXPECT_TRUE(data.ID_passive().isApprox(IDC_ref.head(dim_passive)));
  EXPECT_TRUE(data.ID().isApprox(IDC_ref.segment(dim_passive, dimu)));
  EXPECT_TRUE(data.C().isApprox(IDC_ref.segment(dimv, dimf)));
  EXPECT_TRUE(data.MJtJinv_IDC().isApprox(MJtJinv_IDC_ref));
  EXPECT_TRUE(data.laf().isApprox(laf_ref));
  EXPECT_TRUE(data.la().isApprox(laf_ref.head(dimv)));
  EXPECT_TRUE(data.lf().isApprox(laf_ref.tail(dimf)));
}


TEST_F(ContactDynamicsDataTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(fixed_base_urdf, contact_frames);
  contact_status.setContactStatus({false});
  testSize(robot, contact_status);
  contact_status.setContactStatus({true});
  testSize(robot, contact_status);
}


TEST_F(ContactDynamicsDataTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  ContactStatus contact_status(contact_frames.size());
  Robot robot(floating_base_urdf, contact_frames);
  contact_status.setContactStatus({false, false, false, false});
  testSize(robot, contact_status);
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  testSize(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}