#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_matrix_inverter.hpp"


namespace idocp {

class RiccatiMatrixInverterTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(RiccatiMatrixInverterTest, fixed_base_without_contacts) {
  const int dimv = fixed_base_robot_.dimv();
  const int dimaf = fixed_base_robot_.dimv() + fixed_base_robot_.dimf();
  const int dimf = fixed_base_robot_.dimf();
  const int dimc = fixed_base_robot_.dim_passive() + fixed_base_robot_.dimf();
  ASSERT_EQ(dimf, 0);
  ASSERT_EQ(dimc, 0);
  const Eigen::MatrixXd Qaa_seed = Eigen::MatrixXd::Random(dimaf, dimaf);
  const Eigen::MatrixXd Qaa = Qaa_seed * Qaa_seed.transpose() + Eigen::MatrixXd::Identity(dimaf, dimaf);
  const Eigen::MatrixXd Caf = Eigen::MatrixXd::Random(dimc, dimaf);
  RiccatiMatrixInverter inverter(fixed_base_robot_);
  inverter.setContactStatus(fixed_base_robot_);
  inverter.invert(Qaa, Caf);
  Eigen::MatrixXd G_inv = Eigen::MatrixXd::Zero(dimaf, dimaf);
  inverter.getInverseMatrix(G_inv);
  Eigen::MatrixXd G_inv_ref = Qaa.inverse();
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  EXPECT_TRUE((G_inv*Qaa).isApprox(Eigen::MatrixXd::Identity(dimaf+dimc, dimaf+dimc), 1.0e-12));
  const Eigen::MatrixXd dPvv = Eigen::MatrixXd::Random(dimv, dimv);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  Eigen::MatrixXd dP_ref = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  dP_ref.topLeftCorner(dimv, dimv) = dtau * dtau * dPvv;
  G_inv_ref += G_inv_ref * dP_ref * G_inv_ref;
  inverter.firstOrderCorrection(dtau, dPvv, G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  std::cout << "error l2 norm = " << (G_inv_ref - G_inv).lpNorm<2>() << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, fixed_base_with_contacts) {
  std::vector<int> contact_frames = {18};
  fixed_base_robot_ = Robot(fixed_base_urdf_, contact_frames, 0, 0);
  std::vector<bool> contact_status = {true};
  fixed_base_robot_.setContactStatus(contact_status);
  const int dimv = fixed_base_robot_.dimv();
  const int dimaf = fixed_base_robot_.dimv() + fixed_base_robot_.dimf();
  const int dimf = fixed_base_robot_.dimf();
  const int dimc = fixed_base_robot_.dim_passive() + fixed_base_robot_.dimf();
  const Eigen::MatrixXd Qaa_seed = Eigen::MatrixXd::Random(dimaf, dimaf);
  const Eigen::MatrixXd Qaa = Qaa_seed * Qaa_seed.transpose() + Eigen::MatrixXd::Identity(dimaf, dimaf);
  const Eigen::MatrixXd Caf = Eigen::MatrixXd::Random(dimc, dimaf);
  RiccatiMatrixInverter inverter(fixed_base_robot_);
  inverter.setContactStatus(fixed_base_robot_);
  inverter.invert(Qaa, Caf);
  Eigen::MatrixXd G_mat = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  G_mat.topLeftCorner(dimaf, dimaf) = Qaa;
  G_mat.topRightCorner(dimaf, dimc) = Caf.transpose();
  G_mat.bottomLeftCorner(dimc, dimaf) = Caf;
  std::cout << G_mat << std::endl;
  Eigen::MatrixXd G_inv_ref = G_mat.inverse();
  Eigen::MatrixXd G_inv = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  inverter.getInverseMatrix(G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  EXPECT_TRUE((G_inv*G_mat).isApprox(Eigen::MatrixXd::Identity(dimaf+dimc, dimaf+dimc), 1.0e-12));
  const Eigen::MatrixXd dPvv = Eigen::MatrixXd::Random(dimv, dimv);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  Eigen::MatrixXd dP_ref = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  dP_ref.topLeftCorner(dimv, dimv) = dtau * dtau * dPvv;
  G_inv_ref += G_inv_ref * dP_ref * G_inv_ref;
  inverter.firstOrderCorrection(dtau, dPvv, G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  std::cout << "error l2 norm = " << (G_inv_ref - G_inv).lpNorm<2>() << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, floating_base_without_contacts) {
  const int dimv = floating_base_robot_.dimv();
  const int dimaf = floating_base_robot_.dimv() + floating_base_robot_.dimf();
  const int dimf = floating_base_robot_.dimf();
  const int dimc = floating_base_robot_.dim_passive() + floating_base_robot_.dimf();
  const Eigen::MatrixXd Qaa_seed = Eigen::MatrixXd::Random(dimaf, dimaf);
  const Eigen::MatrixXd Qaa = Qaa_seed * Qaa_seed.transpose() + Eigen::MatrixXd::Identity(dimaf, dimaf);
  const Eigen::MatrixXd Caf = Eigen::MatrixXd::Random(dimc, dimaf);
  RiccatiMatrixInverter inverter(floating_base_robot_);
  inverter.setContactStatus(floating_base_robot_);
  inverter.invert(Qaa, Caf);
  Eigen::MatrixXd G_mat = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  G_mat.topLeftCorner(dimaf, dimaf) = Qaa;
  G_mat.topRightCorner(dimaf, dimc) = Caf.transpose();
  G_mat.bottomLeftCorner(dimc, dimaf) = Caf;
  std::cout << G_mat << std::endl;
  Eigen::MatrixXd G_inv_ref = G_mat.inverse();
  Eigen::MatrixXd G_inv = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  inverter.getInverseMatrix(G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  EXPECT_TRUE((G_inv*G_mat).isApprox(Eigen::MatrixXd::Identity(dimaf+dimc, dimaf+dimc), 1.0e-12));
  std::cout << G_inv_ref - G_inv << std::endl;
  const Eigen::MatrixXd dPvv = Eigen::MatrixXd::Random(dimv, dimv);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  Eigen::MatrixXd dP_ref = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  dP_ref.topLeftCorner(dimv, dimv) = dtau * dtau * dPvv;
  G_inv_ref += G_inv_ref * dP_ref * G_inv_ref;
  inverter.firstOrderCorrection(dtau, dPvv, G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  std::cout << "error l2 norm = " << (G_inv_ref - G_inv).lpNorm<2>() << std::endl;
}


TEST_F(RiccatiMatrixInverterTest, floating_base_with_contacts) {
  const std::vector<int> contact_frames = {14, 24, 34, 44};
  floating_base_robot_ = Robot(floating_base_urdf_, contact_frames, 0, 0);
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<contact_frames.size(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  floating_base_robot_.setContactStatus(active_contacts);
  const int dimv = floating_base_robot_.dimv();
  const int dimaf = floating_base_robot_.dimv() + floating_base_robot_.dimf();
  const int dimf = floating_base_robot_.dimf();
  const int dimc = floating_base_robot_.dim_passive() + floating_base_robot_.dimf();
  const Eigen::MatrixXd Qaa_seed = Eigen::MatrixXd::Random(dimaf, dimaf);
  const Eigen::MatrixXd Qaa = Qaa_seed * Qaa_seed.transpose() + Eigen::MatrixXd::Identity(dimaf, dimaf);
  const Eigen::MatrixXd Caf = Eigen::MatrixXd::Random(dimc, dimaf);
  RiccatiMatrixInverter inverter(floating_base_robot_);
  inverter.setContactStatus(floating_base_robot_);
  inverter.invert(Qaa, Caf);
  Eigen::MatrixXd G_mat = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  G_mat.topLeftCorner(dimaf, dimaf) = Qaa;
  G_mat.topRightCorner(dimaf, dimc) = Caf.transpose();
  G_mat.bottomLeftCorner(dimc, dimaf) = Caf;
  std::cout << G_mat << std::endl;
  Eigen::MatrixXd G_inv_ref = G_mat.inverse();
  Eigen::MatrixXd G_inv = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  inverter.getInverseMatrix(G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  EXPECT_TRUE((G_inv*G_mat).isApprox(Eigen::MatrixXd::Identity(dimaf+dimc, dimaf+dimc), 1.0e-12));
  std::cout << G_inv_ref - G_inv << std::endl;
  const Eigen::MatrixXd dPvv = Eigen::MatrixXd::Random(dimv, dimv);
  const double dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  Eigen::MatrixXd dP_ref = Eigen::MatrixXd::Zero(dimaf+dimc, dimaf+dimc);
  dP_ref.topLeftCorner(dimv, dimv) = dtau * dtau * dPvv;
  G_inv_ref += G_inv_ref * dP_ref * G_inv_ref;
  inverter.firstOrderCorrection(dtau, dPvv, G_inv);
  EXPECT_TRUE(G_inv_ref.isApprox(G_inv, 1.0e-12));
  std::cout << "error l2 norm = " << (G_inv_ref - G_inv).lpNorm<2>() << std::endl;
  std::cout << "dimf = " << dimf << std::endl;
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}