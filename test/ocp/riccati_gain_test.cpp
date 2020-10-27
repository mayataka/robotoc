#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"


namespace idocp {

class RiccatiGainTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
  }

  virtual void TearDown() {
  }

  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(RiccatiGainTest, fixed_base) {
  const int dimv = fixed_base_robot_.dimv();
  const int dimu = fixed_base_robot_.dimu();
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  const Eigen::MatrixXd G_seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd G = G_seed * G_seed.transpose() + Eigen::MatrixXd::Identity(dimu, dimu);
  const Eigen::MatrixXd Qxu = Eigen::MatrixXd::Random(2*dimv, dimu);
  const Eigen::VectorXd lu = Eigen::VectorXd::Random(dimu);
  kkt_matrix.Quu() = G;
  kkt_matrix.Qxu() = Qxu;
  kkt_residual.lu() = lu;
  RiccatiGain gain(fixed_base_robot_);
  gain.computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual);
  const Eigen::MatrixXd Ginv = G.inverse();
  const Eigen::MatrixXd K_ref = - Ginv * Qxu.transpose();
  const Eigen::VectorXd k_ref = - Ginv * lu;
  EXPECT_TRUE(Ginv.isApprox(gain.Ginv));
  EXPECT_TRUE(K_ref.isApprox(gain.K));
  EXPECT_TRUE(k_ref.isApprox(gain.k));
  EXPECT_TRUE(gain.K.leftCols(dimv).isApprox(gain.Kq()));
  EXPECT_TRUE(gain.K.rightCols(dimv).isApprox(gain.Kv()));
}


TEST_F(RiccatiGainTest, floating_base) {
  const int dimv = floating_base_robot_.dimv();
  const int dimu = floating_base_robot_.dimu();
  KKTMatrix kkt_matrix(floating_base_robot_);
  KKTResidual kkt_residual(floating_base_robot_);
  const Eigen::MatrixXd G_seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd G = G_seed * G_seed.transpose() + Eigen::MatrixXd::Identity(dimu, dimu);
  const Eigen::MatrixXd Qxu = Eigen::MatrixXd::Random(2*dimv, dimu);
  const Eigen::VectorXd lu = Eigen::VectorXd::Random(dimu);
  kkt_matrix.Quu() = G;
  kkt_matrix.Qxu() = Qxu;
  kkt_residual.lu() = lu;
  RiccatiGain gain(floating_base_robot_);
  gain.computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual);
  const Eigen::MatrixXd Ginv = G.inverse();
  const Eigen::MatrixXd K_ref = - Ginv * Qxu.transpose();
  const Eigen::VectorXd k_ref = - Ginv * lu;
  EXPECT_TRUE(Ginv.isApprox(gain.Ginv));
  EXPECT_TRUE(K_ref.isApprox(gain.K));
  EXPECT_TRUE(k_ref.isApprox(gain.k));
  EXPECT_TRUE(gain.K.leftCols(dimv).isApprox(gain.Kq()));
  EXPECT_TRUE(gain.K.rightCols(dimv).isApprox(gain.Kv()));
}



} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}