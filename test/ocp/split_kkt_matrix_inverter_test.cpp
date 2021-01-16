#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_matrix_inverter.hpp"


namespace idocp {

class SplitKKTMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void test(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
};


void SplitKKTMatrixInverterTest::test(const Robot& robot) {
  SplitKKTMatrix matrix(robot);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimx = 2*robot.dimv();
  const int dimQ = 2*robot.dimv() + robot.dimu();
  const int dimKKT = 4*robot.dimv() + robot.dimu();
  const Eigen::MatrixXd KKT_seed_mat = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd KKT_mat = KKT_seed_mat * KKT_seed_mat.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  KKT_mat.topLeftCorner(dimx, dimx).setZero();
  matrix.Fqu() = KKT_mat.block(             0,           dimx, dimv, dimu);
  matrix.Fqq() = KKT_mat.block(             0,      dimx+dimu, dimv, dimv);
  matrix.Fqv() = KKT_mat.block(             0, dimx+dimu+dimv, dimv, dimv);
  matrix.Fvu() = KKT_mat.block(          dimv,           dimx, dimv, dimu);
  matrix.Fvq() = KKT_mat.block(          dimv,      dimx+dimu, dimv, dimv);
  matrix.Fvv() = KKT_mat.block(          dimv, dimx+dimu+dimv, dimv, dimv);
  matrix.Quu() = KKT_mat.block(          dimx,           dimx, dimu, dimu);
  matrix.Quq() = KKT_mat.block(          dimx,      dimx+dimu, dimu, dimv);
  matrix.Quv() = KKT_mat.block(          dimx, dimx+dimu+dimv, dimu, dimv);
  matrix.Qqu() = KKT_mat.block(     dimx+dimu,           dimx, dimv, dimu);
  matrix.Qqq() = KKT_mat.block(     dimx+dimu,      dimx+dimu, dimv, dimv);
  matrix.Qqv() = KKT_mat.block(     dimx+dimu, dimx+dimu+dimv, dimv, dimv);
  matrix.Qvu() = KKT_mat.block(dimx+dimu+dimv,           dimx, dimv, dimu);
  matrix.Qvq() = KKT_mat.block(dimx+dimu+dimv,      dimx+dimu, dimv, dimv);
  matrix.Qvv() = KKT_mat.block(dimx+dimu+dimv, dimx+dimu+dimv, dimv, dimv);
  matrix.symmetrize();
  Eigen::MatrixXd KKT_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  SplitKKTMatrixInverter inverter(robot);
  inverter.invert(matrix.Jac(), matrix.Qss(), KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity());
}


TEST_F(SplitKKTMatrixInverterTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus({false});
  test(robot);
  contact_status.setContactStatus({true});
  test(robot);
}


TEST_F(SplitKKTMatrixInverterTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ContactStatus contact_status(contact_frames.size());
  contact_status.setContactStatus(is_contact_active);
  test(robot);
  is_contact_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  test(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}