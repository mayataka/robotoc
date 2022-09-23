#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class ContactDynamicsDataTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }

  static void testSize(const Robot& robot, const ContactStatus& contact_status);

};


void ContactDynamicsDataTest::testSize(const Robot& robot, 
                                       const ContactStatus& contact_status) {
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dim_passive = robot.dim_passive();
  ContactDynamicsData data(robot);
  EXPECT_EQ(data.hasFloatingBase(), robot.hasFloatingBase());
  EXPECT_EQ(data.dimv(), robot.dimv());
  EXPECT_EQ(data.dimvf(), robot.dimv());
  EXPECT_EQ(data.dimu(), robot.dimu());
  data.setContactDimension(contact_status.dimf());
  EXPECT_EQ(data.dimf(), contact_status.dimf());
  EXPECT_EQ(data.dimvf(), robot.dimv()+contact_status.dimf());
  EXPECT_EQ(data.Qxu_passive.rows(), dimx);
  EXPECT_EQ(data.Qxu_passive.cols(), dim_passive);
  EXPECT_EQ(data.Quu_passive_topRight.rows(), dim_passive);
  EXPECT_EQ(data.Quu_passive_topRight.cols(), dimu);
  EXPECT_EQ(data.lu_passive.size(), dim_passive);
  EXPECT_EQ(data.dIDda.rows(), dimv);
  EXPECT_EQ(data.dIDda.cols(), dimv);
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
  EXPECT_EQ(data.haf().size(), dimv+dimf);
  EXPECT_EQ(data.ha().size(), dimv);
  EXPECT_EQ(data.hf().size(), dimf);
  const Eigen::MatrixXd dIDda_ref = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd dCda_ref = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd dIDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd MJtJinv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimv+dimf);
  const Eigen::MatrixXd MJtJinv_dIDCdqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd Qafqv_ref = Eigen::MatrixXd::Random(dimv+dimf, dimx);
  const Eigen::MatrixXd Qafu_full_ref = Eigen::MatrixXd::Random(dimv+dimf, dimv);
  const Eigen::VectorXd IDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd MJtJinv_IDC_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd laf_ref = Eigen::VectorXd::Random(dimv+dimf);
  const Eigen::VectorXd haf_ref = Eigen::VectorXd::Random(dimv+dimf);
  data.dIDda = dIDda_ref;
  data.dCda() = dCda_ref;
  data.dIDCdqv() = dIDCdqv_ref;
  data.MJtJinv() = MJtJinv_ref;
  data.MJtJinv_dIDCdqv() = MJtJinv_dIDCdqv_ref;
  data.Qafqv() = Qafqv_ref;
  data.Qafu_full() = Qafu_full_ref;
  data.IDC() = IDC_ref;
  data.MJtJinv_IDC() = MJtJinv_IDC_ref;
  data.laf() = laf_ref;
  data.haf() = haf_ref;
  EXPECT_TRUE(data.dIDda.isApprox(dIDda_ref));
  EXPECT_TRUE(data.dCda().isApprox(dCda_ref));
  EXPECT_TRUE(data.dIDCdqv().isApprox(dIDCdqv_ref));
  EXPECT_TRUE(data.dIDdq().isApprox(dIDCdqv_ref.topLeftCorner(dimv, dimv)));
  EXPECT_TRUE(data.dIDdv().isApprox(dIDCdqv_ref.topRightCorner(dimv, dimv)));
  EXPECT_TRUE(data.dCdq().isApprox(dIDCdqv_ref.bottomLeftCorner(dimf, dimv)));
  EXPECT_TRUE(data.dCdv().isApprox(dIDCdqv_ref.bottomRightCorner(dimf, dimv)));
  EXPECT_TRUE(data.MJtJinv().isApprox(MJtJinv_ref));
  EXPECT_TRUE(data.MJtJinv_dIDCdqv().isApprox(MJtJinv_dIDCdqv_ref));
  EXPECT_TRUE(data.Qafqv().isApprox(Qafqv_ref));
  EXPECT_TRUE(data.Qafqv().isApprox(data.Qdvfqv()));
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
  EXPECT_TRUE(data.laf().isApprox(data.ldvf()));
  EXPECT_TRUE(data.la().isApprox(data.ldv()));
  EXPECT_TRUE(data.haf().isApprox(haf_ref));
  EXPECT_TRUE(data.ha().isApprox(haf_ref.head(dimv)));
  EXPECT_TRUE(data.hf().isApprox(haf_ref.tail(dimf)));
}


TEST_F(ContactDynamicsDataTest, fixedBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateRobotManipulator(dt);
  auto contact_status = robot.createContactStatus();
  testSize(robot, contact_status);
  contact_status.activateContact(0);
  testSize(robot, contact_status);
}


TEST_F(ContactDynamicsDataTest, floatingBase) {
  const double dt = 0.01;
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  auto contact_status = robot.createContactStatus();
  testSize(robot, contact_status);
  contact_status.setRandom();
  testSize(robot, contact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}