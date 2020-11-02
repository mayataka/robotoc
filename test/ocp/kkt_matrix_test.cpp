#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class KKTMatrixTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
  }

  virtual void TearDown() {
  }

  static void testSize(const Robot& robot, const ContactStatus& contact_status, const ImpulseStatus& impulse_status);
  static void testIsApprox(const Robot& robot, const ContactStatus& contact_status, const ImpulseStatus& impulse_status);
  static void testInverse(const Robot& robot, const ImpulseStatus& impulse_status);

  double dtau_;
  std::string fixed_base_urdf, floating_base_urdf;
};


void KKTMatrixTest::testSize(const Robot& robot, const ContactStatus& contact_status, const ImpulseStatus& impulse_status) {
  KKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dim_passive = robot.dim_passive();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimKKT(), 4*dimv+dimu);
  EXPECT_EQ(matrix.Fqu().rows(), dimv);
  EXPECT_EQ(matrix.Fqu().cols(), dimu);
  EXPECT_EQ(matrix.Fqq().rows(), dimv);
  EXPECT_EQ(matrix.Fqq().cols(), dimv);
  EXPECT_EQ(matrix.Fqv().rows(), dimv);
  EXPECT_EQ(matrix.Fqv().cols(), dimv);
  EXPECT_EQ(matrix.Fvu().rows(), dimv);
  EXPECT_EQ(matrix.Fvu().cols(), dimu);
  EXPECT_EQ(matrix.Fvq().rows(), dimv);
  EXPECT_EQ(matrix.Fvq().cols(), dimv);
  EXPECT_EQ(matrix.Fvv().rows(), dimv);
  EXPECT_EQ(matrix.Fvv().cols(), dimv);
  EXPECT_EQ(matrix.Fxu().rows(), 2*dimv);
  EXPECT_EQ(matrix.Fxu().cols(), dimu);
  EXPECT_EQ(matrix.Fxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Fxx().cols(), 2*dimv);
  EXPECT_EQ(matrix.Quu_full().rows(), dimv);
  EXPECT_EQ(matrix.Quu_full().cols(), dimv);
  EXPECT_EQ(matrix.Quu_passive_topLeft().rows(), dim_passive);
  EXPECT_EQ(matrix.Quu_passive_topLeft().cols(), dim_passive);
  EXPECT_EQ(matrix.Quu_passive_topRight().rows(), dim_passive);
  EXPECT_EQ(matrix.Quu_passive_topRight().cols(), dimu);
  EXPECT_EQ(matrix.Quu_passive_bottomLeft().rows(), dimu);
  EXPECT_EQ(matrix.Quu_passive_bottomLeft().cols(), dim_passive);
  EXPECT_EQ(matrix.Quu().rows(), dimu);
  EXPECT_EQ(matrix.Quu().cols(), dimu);
  EXPECT_EQ(matrix.Quq_full().rows(), dimv);
  EXPECT_EQ(matrix.Quq_full().cols(), dimv);
  EXPECT_EQ(matrix.Quq_passive().rows(), dim_passive);
  EXPECT_EQ(matrix.Quq_passive().cols(), dimv);
  EXPECT_EQ(matrix.Quq().rows(), dimu);
  EXPECT_EQ(matrix.Quq().cols(), dimv);
  EXPECT_EQ(matrix.Quv_full().rows(), dimv);
  EXPECT_EQ(matrix.Quv_full().cols(), dimv);
  EXPECT_EQ(matrix.Quv_passive().rows(), dim_passive);
  EXPECT_EQ(matrix.Quv_passive().cols(), dimv);
  EXPECT_EQ(matrix.Quv().rows(), dimu);
  EXPECT_EQ(matrix.Quv().cols(), dimv);
  EXPECT_EQ(matrix.Qqu_full().rows(), dimv);
  EXPECT_EQ(matrix.Qqu_full().cols(), dimv);
  EXPECT_EQ(matrix.Qqu_passive().rows(), dimv);
  EXPECT_EQ(matrix.Qqu_passive().cols(), dim_passive);
  EXPECT_EQ(matrix.Qqu().rows(), dimv);
  EXPECT_EQ(matrix.Qqu().cols(), dimu);
  EXPECT_EQ(matrix.Qqq().rows(), dimv);
  EXPECT_EQ(matrix.Qqq().cols(), dimv);
  EXPECT_EQ(matrix.Qqv().rows(), dimv);
  EXPECT_EQ(matrix.Qqv().cols(), dimv);
  EXPECT_EQ(matrix.Qvu_full().rows(), dimv);
  EXPECT_EQ(matrix.Qvu_full().cols(), dimv);
  EXPECT_EQ(matrix.Qvu_passive().rows(), dimv);
  EXPECT_EQ(matrix.Qvu_passive().cols(), dim_passive);
  EXPECT_EQ(matrix.Qvu().rows(), dimv);
  EXPECT_EQ(matrix.Qvu().cols(), dimu);
  EXPECT_EQ(matrix.Qvq().rows(), dimv);
  EXPECT_EQ(matrix.Qvq().cols(), dimv);
  EXPECT_EQ(matrix.Qvv().rows(), dimv);
  EXPECT_EQ(matrix.Qvv().cols(), dimv);
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_EQ(matrix.Qxu_full().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxu_full().cols(), dimv);
  EXPECT_EQ(matrix.Qxu_passive().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxu_passive().cols(), dim_passive);
  EXPECT_EQ(matrix.Qxu().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxu().cols(), dimu);
  EXPECT_EQ(matrix.Qux_full().rows(), dimv);
  EXPECT_EQ(matrix.Qux_full().cols(), 2*dimv);
  EXPECT_EQ(matrix.Qux_passive().rows(), dim_passive);
  EXPECT_EQ(matrix.Qux_passive().cols(), 2*dimv);
  EXPECT_EQ(matrix.Qux().rows(), dimu);
  EXPECT_EQ(matrix.Qux().cols(), 2*dimv);
  EXPECT_EQ(matrix.Qaa().rows(), dimv);
  EXPECT_EQ(matrix.Qaa().cols(), dimv);
  EXPECT_EQ(matrix.Qff().rows(), dimf);
  EXPECT_EQ(matrix.Qff().cols(), dimf);
  EXPECT_EQ(matrix.Qaaff().rows(), dimv+dimf);
  EXPECT_EQ(matrix.Qaaff().cols(), dimv+dimf);
  EXPECT_EQ(matrix.Fqq_prev.rows(), dimv);
  EXPECT_EQ(matrix.Fqq_prev.cols(), dimv);
  const Eigen::MatrixXd Fqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quq_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quv_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Fqu() = Fqu;
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvu() = Fvu;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
  matrix.Quu_full() = Quu_full;
  matrix.Quq_full() = Quq_full;
  matrix.Quv_full() = Quv_full;
  matrix.Qqu_full() = Qqu_full;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvu_full() = Qvu_full;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  matrix.Qff() = Qff;
  matrix.Qaa() = Qaa;
  EXPECT_TRUE(matrix.Fqu().isApprox(Fqu));
  EXPECT_TRUE(matrix.Fqq().isApprox(Fqq));
  EXPECT_TRUE(matrix.Fqv().isApprox(Fqv));
  EXPECT_TRUE(matrix.Fvu().isApprox(Fvu));
  EXPECT_TRUE(matrix.Fvq().isApprox(Fvq));
  EXPECT_TRUE(matrix.Fvv().isApprox(Fvv));
  EXPECT_TRUE(matrix.Fxu().topRows(dimv).isApprox(Fqu));
  EXPECT_TRUE(matrix.Fxu().bottomRows(dimv).isApprox(Fvu));
  EXPECT_TRUE(matrix.Fxx().topLeftCorner(dimv, dimv).isApprox(Fqq));
  EXPECT_TRUE(matrix.Fxx().topRightCorner(dimv, dimv).isApprox(Fqv));
  EXPECT_TRUE(matrix.Fxx().bottomLeftCorner(dimv, dimv).isApprox(Fvq));
  EXPECT_TRUE(matrix.Fxx().bottomRightCorner(dimv, dimv).isApprox(Fvv));
  EXPECT_TRUE(matrix.Quu_full().isApprox(Quu_full));
  EXPECT_TRUE(matrix.Quu().isApprox(Quu_full.bottomRightCorner(dimu, dimu)));
  EXPECT_TRUE(matrix.Quu_passive_topLeft().isApprox(Quu_full.topLeftCorner(dim_passive, dim_passive)));
  EXPECT_TRUE(matrix.Quu_passive_topRight().isApprox(Quu_full.topRightCorner(dim_passive, dimu)));
  EXPECT_TRUE(matrix.Quu_passive_bottomLeft().isApprox(Quu_full.bottomLeftCorner(dimu, dim_passive)));
  EXPECT_TRUE(matrix.Quq_full().isApprox(Quq_full));
  EXPECT_TRUE(matrix.Quq_passive().isApprox(Quq_full.topRows(dim_passive)));
  EXPECT_TRUE(matrix.Quq().isApprox(Quq_full.bottomRows(dimu)));
  EXPECT_TRUE(matrix.Quv_full().isApprox(Quv_full));
  EXPECT_TRUE(matrix.Quv_passive().isApprox(Quv_full.topRows(dim_passive)));
  EXPECT_TRUE(matrix.Quv().isApprox(Quv_full.bottomRows(dimu)));
  EXPECT_TRUE(matrix.Qqu_full().isApprox(Qqu_full));
  EXPECT_TRUE(matrix.Qqu_passive().isApprox(Qqu_full.leftCols(dim_passive)));
  EXPECT_TRUE(matrix.Qqu().isApprox(Qqu_full.rightCols(dimu)));
  EXPECT_TRUE(matrix.Qqq().isApprox(Qqq));
  EXPECT_TRUE(matrix.Qqv().isApprox(Qqv));
  EXPECT_TRUE(matrix.Qvu_full().isApprox(Qvu_full));
  EXPECT_TRUE(matrix.Qvu_passive().isApprox(Qvu_full.leftCols(dim_passive)));
  EXPECT_TRUE(matrix.Qvu().isApprox(Qvu_full.rightCols(dimu)));
  EXPECT_TRUE(matrix.Qvq().isApprox(Qvq));
  EXPECT_TRUE(matrix.Qvv().isApprox(Qvv));
  EXPECT_TRUE(matrix.Qaa().isApprox(Qaa));
  EXPECT_TRUE(matrix.Qff().isApprox(Qff));
  EXPECT_TRUE(matrix.Qaa().isApprox(Qaa));
  EXPECT_TRUE(matrix.Qff().isApprox(Qff));
  EXPECT_TRUE(matrix.Qaaff().topLeftCorner(dimv, dimv).isApprox(Qaa));
  EXPECT_TRUE(matrix.Qaaff().bottomRightCorner(dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.Qxx().topLeftCorner(dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().topRightCorner(dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().bottomLeftCorner(dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().bottomRightCorner(dimv, dimv).isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxu_full().topRows(dimv).isApprox(Qqu_full));
  EXPECT_TRUE(matrix.Qxu_full().bottomRows(dimv).isApprox(Qvu_full));
  EXPECT_TRUE(matrix.Qxu_passive().topRows(dimv).isApprox(Qqu_full.leftCols(dim_passive)));
  EXPECT_TRUE(matrix.Qxu_passive().bottomRows(dimv).isApprox(Qvu_full.leftCols(dim_passive)));
  EXPECT_TRUE(matrix.Qxu().topRows(dimv).isApprox(Qqu_full.rightCols(dimu)));
  EXPECT_TRUE(matrix.Qxu().bottomRows(dimv).isApprox(Qvu_full.rightCols(dimu)));
  EXPECT_TRUE(matrix.Qxu_passive().isApprox(matrix.Qxu_full().leftCols(dim_passive)));
  EXPECT_TRUE(matrix.Qxu().isApprox(matrix.Qxu_full().rightCols(dimu)));
  EXPECT_TRUE(matrix.Qux_full().leftCols(dimv).isApprox(Quq_full));
  EXPECT_TRUE(matrix.Qux_full().rightCols(dimv).isApprox(Quv_full));
  EXPECT_TRUE(matrix.Qux_passive().leftCols(dimv).isApprox(Quq_full.topRows(dim_passive)));
  EXPECT_TRUE(matrix.Qux_passive().rightCols(dimv).isApprox(Quv_full.topRows(dim_passive)));
  EXPECT_TRUE(matrix.Qux().leftCols(dimv).isApprox(Quq_full.bottomRows(dimu)));
  EXPECT_TRUE(matrix.Qux().rightCols(dimv).isApprox(Quv_full.bottomRows(dimu)));
  EXPECT_TRUE(matrix.Qux_passive().isApprox(matrix.Qux_full().topRows(dim_passive)));
  EXPECT_TRUE(matrix.Qux().isApprox(matrix.Qux_full().bottomRows(dimu)));
}


void KKTMatrixTest::testIsApprox(const Robot& robot, const ContactStatus& contact_status, const ImpulseStatus& impulse_status) {
  KKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const int dimf = contact_status.dimf();
  const int dim_passive = robot.dim_passive();
  const Eigen::MatrixXd Fqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Fqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Fvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Fvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quq_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Quv_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvu_full = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qaa = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Fqu() = Fqu;
  matrix.Fqq() = Fqq;
  matrix.Fqv() = Fqv;
  matrix.Fvu() = Fvu;
  matrix.Fvq() = Fvq;
  matrix.Fvv() = Fvv;
  matrix.Quu_full() = Quu_full;
  matrix.Quq_full() = Quq_full;
  matrix.Quv_full() = Quv_full;
  matrix.Qqu_full() = Qqu_full;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvu_full() = Qvu_full;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  matrix.Qff() = Qff;
  matrix.Qaa() = Qaa;
  KKTMatrix matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Fxx().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Fxu().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Qxx().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Qxu().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Qux().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
  matrix_ref.Qaaff().setRandom();
  EXPECT_FALSE(matrix.isApprox(matrix_ref));
  matrix_ref = matrix;
  EXPECT_TRUE(matrix.isApprox(matrix_ref));
}


void KKTMatrixTest::testInverse(const Robot& robot, const ImpulseStatus& impulse_status) {
  KKTMatrix matrix(robot);
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
  matrix.invert(KKT_mat_inv);
  const Eigen::MatrixXd KKT_mat_inv_ref = KKT_mat.inverse();
  EXPECT_TRUE(KKT_mat_inv.isApprox(KKT_mat_inv_ref));
  EXPECT_TRUE((KKT_mat_inv*KKT_mat).isIdentity());
}


TEST_F(KKTMatrixTest, fixedBase) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  ImpulseStatus impulse_status(contact_frames.size());
  contact_status.setContactStatus({false});
  impulse_status.setImpulseStatus({false});
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  contact_status.setContactStatus({true});
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  testInverse(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  testInverse(robot, impulse_status);
}


TEST_F(KKTMatrixTest, floatingBase) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, contact_frames);
  std::vector<bool> is_contact_active = {false, false, false, false};
  ContactStatus contact_status(contact_frames.size());
  ImpulseStatus impulse_status(contact_frames.size());
  contact_status.setContactStatus(is_contact_active);
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  is_contact_active.clear();
  std::random_device rnd;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  testInverse(robot, impulse_status);
  impulse_status.setImpulseStatus(is_contact_active);
  testSize(robot, contact_status, impulse_status);
  testIsApprox(robot, contact_status, impulse_status);
  testInverse(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}