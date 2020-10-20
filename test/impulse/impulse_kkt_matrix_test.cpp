#include <string>

#include <gtest/gtest.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

class ImpulseKKTMatrixTest : public ::testing::Test {
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


TEST_F(ImpulseKKTMatrixTest, fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  ImpulseKKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimc(), dimc);
  EXPECT_EQ(matrix.dimKKT(), 4*dimv+dimf+dimc);
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(dimc, dimf);
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Cf() = Cf;
  matrix.Cq() = Cq;
  matrix.Cv() = Cv;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qfv() = Qfv;
  matrix.Qqf() = Qqf;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvf() = Qvf;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,         0, dimc, dimf).isApprox(Cf));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,      dimf, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0, dimv+dimf, dimc, dimv).isApprox(Cv));
  EXPECT_TRUE(matrix.costHessian().block(        0,           0, dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.costHessian().block(        0,        dimf, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.costHessian().block(        0,   dimf+dimv, dimf, dimv).isApprox(Qfv));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,           0, dimv, dimf).isApprox(Qqf));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,        dimf, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,   dimf+dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,           0, dimv, dimf).isApprox(Qvf));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,        dimf, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,   dimf+dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Cqv().rows(), dimc);
  EXPECT_EQ(matrix.Cqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Cqv().block(0, 0, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.Cqv().block(0, dimv, dimc, dimv).isApprox(Cv));
  matrix.setZero();
  EXPECT_TRUE(matrix.costHessian().isZero());
  EXPECT_TRUE(matrix.constraintsJacobian().isZero());
}


TEST_F(ImpulseKKTMatrixTest, invert_fixed_base) {
  std::vector<int> contact_frames = {18};
  Robot robot(fixed_base_urdf_, contact_frames);
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active = {true};
  contact_status.setContactStatus(is_contact_active);
  ImpulseKKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf();
  const int dimQ = 2*robot.dimv() + contact_status.dimf();
  const int dimKKT = 4*dimv+2*dimf;
  EXPECT_EQ(dimx+dimQ+dimc, dimKKT);
  const Eigen::MatrixXd kkt_mat_seed = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd kkt_mat = kkt_mat_seed * kkt_mat_seed.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  kkt_mat.topLeftCorner(dimx+dimc, dimx+dimc).setZero();
  kkt_mat.block(0, dimx+dimc, dimv, dimf).setZero();
  kkt_mat.block(0, dimx+dimc+dimf, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat.block(0, dimx+dimc+dimf+dimv, dimv, dimv).setZero();
  kkt_mat.block(dimv, dimx+dimc+dimf+dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat.block(dimx, dimx+dimc, dimc, dimf).setZero();
  kkt_mat.block(dimx, dimx+dimc+dimf, dimf, dimv).setZero();
  kkt_mat.triangularView<Eigen::StrictlyLower>() = kkt_mat.transpose().triangularView<Eigen::StrictlyLower>();
  const Eigen::MatrixXd kkt_mat_inv_ref = kkt_mat.inverse();
  matrix.costHessian() = kkt_mat.bottomRightCorner(dimQ, dimQ);
  matrix.constraintsJacobian() = kkt_mat.block(dimx, dimx+dimc, dimc, dimQ);
  matrix.Fqq = kkt_mat.block(0, dimx+dimc+dimf, dimv, dimv);
  matrix.Fvf() = kkt_mat.block(dimv, dimx+dimc, dimv, dimf);
  matrix.Fvq = kkt_mat.block(dimv, dimx+dimc+dimf, dimv, dimv);
  matrix.Fvv = kkt_mat.block(dimv, dimx+dimc+dimf+dimv, dimv, dimv);
  Eigen::MatrixXd kkt_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  matrix.symmetrize();
  matrix.invert(kkt_mat_inv);
  EXPECT_TRUE(kkt_mat_inv.isApprox(kkt_mat_inv_ref, 1.0e-10));
}


TEST_F(ImpulseKKTMatrixTest, floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  std::random_device rnd;
  ContactStatus contact_status(contact_frames.size());
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  ImpulseKKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf();
  EXPECT_EQ(matrix.dimf(), dimf);
  EXPECT_EQ(matrix.dimc(), dimc);
  EXPECT_EQ(matrix.dimKKT(), 4*dimv+dimf+dimc);
  const Eigen::MatrixXd Cf = Eigen::MatrixXd::Random(dimc, dimf);
  const Eigen::MatrixXd Cq = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Cv = Eigen::MatrixXd::Random(dimc, dimv);
  const Eigen::MatrixXd Qff = Eigen::MatrixXd::Random(dimf, dimf);
  const Eigen::MatrixXd Qfq = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qfv = Eigen::MatrixXd::Random(dimf, dimv);
  const Eigen::MatrixXd Qqf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qqq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvf = Eigen::MatrixXd::Random(dimv, dimf);
  const Eigen::MatrixXd Qvq = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = Eigen::MatrixXd::Random(dimv, dimv);
  matrix.Cf() = Cf;
  matrix.Cq() = Cq;
  matrix.Cv() = Cv;
  matrix.Qff() = Qff;
  matrix.Qfq() = Qfq;
  matrix.Qfv() = Qfv;
  matrix.Qqf() = Qqf;
  matrix.Qqq() = Qqq;
  matrix.Qqv() = Qqv;
  matrix.Qvf() = Qvf;
  matrix.Qvq() = Qvq;
  matrix.Qvv() = Qvv;
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,         0, dimc, dimf).isApprox(Cf));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0,      dimf, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.constraintsJacobian().block(0, dimv+dimf, dimc, dimv).isApprox(Cv));
  EXPECT_TRUE(matrix.costHessian().block(        0,           0, dimf, dimf).isApprox(Qff));
  EXPECT_TRUE(matrix.costHessian().block(        0,        dimf, dimf, dimv).isApprox(Qfq));
  EXPECT_TRUE(matrix.costHessian().block(        0,   dimf+dimv, dimf, dimv).isApprox(Qfv));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,           0, dimv, dimf).isApprox(Qqf));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,        dimf, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.costHessian().block(     dimf,   dimf+dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,           0, dimv, dimf).isApprox(Qvf));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,        dimf, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.costHessian().block(dimv+dimf,   dimf+dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Qxx().rows(), 2*dimv);
  EXPECT_EQ(matrix.Qxx().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_TRUE(matrix.Qxx().block(   0,    0, dimv, dimv).isApprox(Qqq));
  EXPECT_TRUE(matrix.Qxx().block(   0, dimv, dimv, dimv).isApprox(Qqv));
  EXPECT_TRUE(matrix.Qxx().block(dimv,    0, dimv, dimv).isApprox(Qvq));
  EXPECT_TRUE(matrix.Qxx().block(dimv, dimv, dimv, dimv).isApprox(Qvv));
  EXPECT_EQ(matrix.Cqv().rows(), dimc);
  EXPECT_EQ(matrix.Cqv().cols(), 2*dimv);
  EXPECT_TRUE(matrix.Cqv().block(0, 0, dimc, dimv).isApprox(Cq));
  EXPECT_TRUE(matrix.Cqv().block(0, dimv, dimc, dimv).isApprox(Cv));
  matrix.setZero();
  EXPECT_TRUE(matrix.costHessian().isZero());
  EXPECT_TRUE(matrix.constraintsJacobian().isZero());
}


TEST_F(ImpulseKKTMatrixTest, invert_floating_base) {
  std::vector<int> contact_frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf_, contact_frames);
  ContactStatus contact_status(contact_frames.size());
  std::random_device rnd;
  std::vector<bool> is_contact_active;
  for (const auto frame : contact_frames) {
    is_contact_active.push_back(rnd()%2==0);
  }
  contact_status.setContactStatus(is_contact_active);
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  ImpulseKKTMatrix matrix(robot);
  matrix.setContactStatus(contact_status);
  const int dimv = robot.dimv();
  const int dimx = 2*robot.dimv();
  const int dimf = contact_status.dimf();
  const int dimc = contact_status.dimf();
  const int dimQ = 2*robot.dimv() + contact_status.dimf();
  const int dimKKT = 4*dimv+2*dimf;
  EXPECT_EQ(dimx+dimQ+dimc, dimKKT);
  const Eigen::MatrixXd kkt_mat_seed = Eigen::MatrixXd::Random(dimKKT, dimKKT);
  Eigen::MatrixXd kkt_mat = kkt_mat_seed * kkt_mat_seed.transpose() + Eigen::MatrixXd::Identity(dimKKT, dimKKT);
  kkt_mat.topLeftCorner(dimx+dimc, dimx+dimc).setZero();
  kkt_mat.block(0, dimx+dimc, dimv, dimf).setZero();
  kkt_mat.block(0, dimx+dimc+dimf, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat.block(0, dimx+dimc+dimf, dimv, dimv).topLeftCorner(6, 6).setRandom();
  kkt_mat.block(0, dimx+dimc+dimf+dimv, dimv, dimv).setZero();
  kkt_mat.block(dimv, dimx+dimc+dimf+dimv, dimv, dimv) = - Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_mat.block(dimx, dimx+dimc, dimc, dimf).setZero();
  kkt_mat.block(dimx, dimx+dimc+dimf, dimf, dimv).setZero();
  kkt_mat.triangularView<Eigen::StrictlyLower>() = kkt_mat.transpose().triangularView<Eigen::StrictlyLower>();
  const Eigen::MatrixXd kkt_mat_inv_ref = kkt_mat.inverse();
  matrix.costHessian() = kkt_mat.bottomRightCorner(dimQ, dimQ);
  matrix.constraintsJacobian() = kkt_mat.block(dimx, dimx+dimc, dimc, dimQ);
  matrix.Fqq = kkt_mat.block(0, dimx+dimc+dimf, dimv, dimv);
  matrix.Fvf() = kkt_mat.block(dimv, dimx+dimc, dimv, dimf);
  matrix.Fvq = kkt_mat.block(dimv, dimx+dimc+dimf, dimv, dimv);
  matrix.Fvv = kkt_mat.block(dimv, dimx+dimc+dimf+dimv, dimv, dimv);
  Eigen::MatrixXd kkt_mat_inv = Eigen::MatrixXd::Zero(dimKKT, dimKKT);
  matrix.symmetrize();
  matrix.invert(kkt_mat_inv);
  EXPECT_TRUE(kkt_mat_inv.isApprox(kkt_mat_inv_ref, 1.0e-08));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}