#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_gain.hpp"


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
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(RiccatiGainTest, fixed_base_without_contacts) {
  const int dimv = fixed_base_robot_.dimv();
  const int dimaf = fixed_base_robot_.dimv() + fixed_base_robot_.dimf();
  const int dimf = fixed_base_robot_.dimf();
  const int dimc = fixed_base_robot_.dim_passive() + fixed_base_robot_.dimf();
  ASSERT_EQ(dimf, 0);
  ASSERT_EQ(dimc, 0);
  const Eigen::MatrixXd Ginv_seed = Eigen::MatrixXd::Random(dimaf+dimc, dimaf+dimc);
  const Eigen::MatrixXd Ginv = Ginv_seed * Ginv_seed.transpose();
  const Eigen::MatrixXd Qafqv = Eigen::MatrixXd::Random(dimaf, 2*dimv);
  const Eigen::MatrixXd Cqv = Eigen::MatrixXd::Random(dimc, 2*dimv);
  RiccatiGain gain(fixed_base_robot_);
  gain.computeFeedbackGain(Ginv, Qafqv, Cqv);
  Eigen::MatrixXd HC = Eigen::MatrixXd::Random(dimaf+dimc, 2*dimv);
  HC.topRows(dimaf) = Qafqv;
  HC.bottomRows(dimc) = Cqv;
  const Eigen::MatrixXd K_ref = - Ginv * HC;
  ASSERT_EQ(K_ref.rows(), dimv);
  ASSERT_EQ(K_ref.cols(), 2*dimv);
  EXPECT_TRUE(K_ref.topLeftCorner(dimv, dimv).isApprox(gain.Kaq()));
  EXPECT_TRUE(K_ref.topRightCorner(dimv, dimv).isApprox(gain.Kav()));
  EXPECT_TRUE(K_ref.block(dimv, 0, dimf, dimv).isApprox(gain.Kfq()));
  EXPECT_TRUE(K_ref.block(dimv, dimv, dimf, dimv).isApprox(gain.Kfv()));
  EXPECT_TRUE(K_ref.bottomLeftCorner(dimc, dimv).isApprox(gain.Kmuq()));
  EXPECT_TRUE(K_ref.bottomRightCorner(dimc, dimv).isApprox(gain.Kmuv()));
  const Eigen::VectorXd laf = Eigen::VectorXd::Random(dimaf);
  const Eigen::VectorXd C = Eigen::VectorXd::Random(dimc);
  gain.computeFeedforward(Ginv, laf, C);
  Eigen::VectorXd h = Eigen::VectorXd::Random(dimaf+dimc);
  h.head(dimaf) = laf;
  h.tail(dimc) = C;
  const Eigen::VectorXd k_ref = - Ginv * h;
  ASSERT_EQ(k_ref.size(), dimv);
  EXPECT_TRUE(k_ref.head(dimv).isApprox(gain.ka()));
  EXPECT_TRUE(k_ref.segment(dimv, dimf).isApprox(gain.kf()));
  EXPECT_TRUE(k_ref.tail(dimc).isApprox(gain.kmu()));
}


TEST_F(RiccatiGainTest, fixed_base_with_contacts) {
  std::vector<int> contact_frames = {18};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  fixed_base_robot_ = Robot(fixed_base_urdf_, contact_frames, mu, 0, 0);
  std::vector<bool> contact_status = {true};
  fixed_base_robot_.setContactStatus(contact_status);
  const int dimv = fixed_base_robot_.dimv();
  const int dimaf = fixed_base_robot_.dimv() + fixed_base_robot_.dimf();
  const int dimf = fixed_base_robot_.dimf();
  const int dimc = fixed_base_robot_.dim_passive() + fixed_base_robot_.dimf();
  ASSERT_EQ(dimf, 3);
  ASSERT_EQ(dimc, 3);
  const Eigen::MatrixXd Ginv_seed = Eigen::MatrixXd::Random(dimaf+dimc, dimaf+dimc);
  const Eigen::MatrixXd Ginv = Ginv_seed * Ginv_seed.transpose();
  const Eigen::MatrixXd Qafqv = Eigen::MatrixXd::Random(dimaf, 2*dimv);
  const Eigen::MatrixXd Cqv = Eigen::MatrixXd::Random(dimc, 2*dimv);
  RiccatiGain gain(fixed_base_robot_);
  gain.computeFeedbackGain(Ginv, Qafqv, Cqv);
  Eigen::MatrixXd HC = Eigen::MatrixXd::Random(dimaf+dimc, 2*dimv);
  HC.topRows(dimaf) = Qafqv;
  HC.bottomRows(dimc) = Cqv;
  const Eigen::MatrixXd K_ref = - Ginv * HC;
  ASSERT_EQ(K_ref.rows(), dimaf+dimc);
  ASSERT_EQ(K_ref.cols(), 2*dimv);
  EXPECT_TRUE(K_ref.topLeftCorner(dimv, dimv).isApprox(gain.Kaq()));
  EXPECT_TRUE(K_ref.topRightCorner(dimv, dimv).isApprox(gain.Kav()));
  EXPECT_TRUE(K_ref.block(dimv, 0, dimf, dimv).isApprox(gain.Kfq()));
  EXPECT_TRUE(K_ref.block(dimv, dimv, dimf, dimv).isApprox(gain.Kfv()));
  EXPECT_TRUE(K_ref.bottomLeftCorner(dimc, dimv).isApprox(gain.Kmuq()));
  EXPECT_TRUE(K_ref.bottomRightCorner(dimc, dimv).isApprox(gain.Kmuv()));
  const Eigen::VectorXd laf = Eigen::VectorXd::Random(dimaf);
  const Eigen::VectorXd C = Eigen::VectorXd::Random(dimc);
  gain.computeFeedforward(Ginv, laf, C);
  Eigen::VectorXd h = Eigen::VectorXd::Random(dimaf+dimc);
  h.head(dimaf) = laf;
  h.tail(dimc) = C;
  const Eigen::VectorXd k_ref = - Ginv * h;
  ASSERT_EQ(k_ref.size(), dimaf+dimc);
  EXPECT_TRUE(k_ref.head(dimv).isApprox(gain.ka()));
  EXPECT_TRUE(k_ref.segment(dimv, dimf).isApprox(gain.kf()));
  EXPECT_TRUE(k_ref.tail(dimc).isApprox(gain.kmu()));
}


TEST_F(RiccatiGainTest, floating_base_without_contacts) {
  const int dimv = floating_base_robot_.dimv();
  const int dimaf = floating_base_robot_.dimv() + floating_base_robot_.dimf();
  const int dimf = floating_base_robot_.dimf();
  const int dimc = floating_base_robot_.dim_passive() + floating_base_robot_.dimf();
  ASSERT_EQ(dimf, 0);
  ASSERT_EQ(dimc, 6);
  const Eigen::MatrixXd Ginv_seed = Eigen::MatrixXd::Random(dimaf+dimc, dimaf+dimc);
  const Eigen::MatrixXd Ginv = Ginv_seed * Ginv_seed.transpose();
  const Eigen::MatrixXd Qafqv = Eigen::MatrixXd::Random(dimaf, 2*dimv);
  const Eigen::MatrixXd Cqv = Eigen::MatrixXd::Random(dimc, 2*dimv);
  RiccatiGain gain(floating_base_robot_);
  gain.computeFeedbackGain(Ginv, Qafqv, Cqv);
  Eigen::MatrixXd HC = Eigen::MatrixXd::Random(dimaf+dimc, 2*dimv);
  HC.topRows(dimaf) = Qafqv;
  HC.bottomRows(dimc) = Cqv;
  const Eigen::MatrixXd K_ref = - Ginv * HC;
  ASSERT_EQ(K_ref.rows(), dimaf+dimc);
  ASSERT_EQ(K_ref.cols(), 2*dimv);
  EXPECT_TRUE(K_ref.topLeftCorner(dimv, dimv).isApprox(gain.Kaq()));
  EXPECT_TRUE(K_ref.topRightCorner(dimv, dimv).isApprox(gain.Kav()));
  EXPECT_TRUE(K_ref.block(dimv, 0, dimf, dimv).isApprox(gain.Kfq()));
  EXPECT_TRUE(K_ref.block(dimv, dimv, dimf, dimv).isApprox(gain.Kfv()));
  EXPECT_TRUE(K_ref.bottomLeftCorner(dimc, dimv).isApprox(gain.Kmuq()));
  EXPECT_TRUE(K_ref.bottomRightCorner(dimc, dimv).isApprox(gain.Kmuv()));
  const Eigen::VectorXd laf = Eigen::VectorXd::Random(dimaf);
  const Eigen::VectorXd C = Eigen::VectorXd::Random(dimc);
  gain.computeFeedforward(Ginv, laf, C);
  Eigen::VectorXd h = Eigen::VectorXd::Random(dimaf+dimc);
  h.head(dimaf) = laf;
  h.tail(dimc) = C;
  const Eigen::VectorXd k_ref = - Ginv * h;
  ASSERT_EQ(k_ref.size(), dimaf+dimc);
  EXPECT_TRUE(k_ref.head(dimv).isApprox(gain.ka()));
  EXPECT_TRUE(k_ref.segment(dimv, dimf).isApprox(gain.kf()));
  EXPECT_TRUE(k_ref.tail(dimc).isApprox(gain.kmu()));
}


TEST_F(RiccatiGainTest, floating_base_with_contacts) {
  const std::vector<int> contact_frames = {14, 24, 34, 44};
  std::vector<double> mu;
  for (int i=0; i<contact_frames.size(); ++i) {
    mu.push_back(std::abs(Eigen::VectorXd::Random(1)[0]));
  }
  floating_base_robot_ = Robot(floating_base_urdf_, contact_frames, mu, 0, 0);
  std::vector<bool> active_contacts;
  std::random_device rnd;
  for (int i=0; i<contact_frames.size(); ++i) {
    active_contacts.push_back(rnd()%1==0);
  }
  floating_base_robot_.setContactStatus(active_contacts);
  const int dimv = floating_base_robot_.dimv();
  const int dimaf = floating_base_robot_.dimv() + floating_base_robot_.dimf();
  const int dimf = floating_base_robot_.dimf();
  const int dimc = floating_base_robot_.dim_passive() + floating_base_robot_.dimf();
  const Eigen::MatrixXd Ginv_seed = Eigen::MatrixXd::Random(dimaf+dimc, dimaf+dimc);
  const Eigen::MatrixXd Ginv = Ginv_seed * Ginv_seed.transpose();
  const Eigen::MatrixXd Qafqv = Eigen::MatrixXd::Random(dimaf, 2*dimv);
  const Eigen::MatrixXd Cqv = Eigen::MatrixXd::Random(dimc, 2*dimv);
  RiccatiGain gain(floating_base_robot_);
  gain.computeFeedbackGain(Ginv, Qafqv, Cqv);
  Eigen::MatrixXd HC = Eigen::MatrixXd::Random(dimaf+dimc, 2*dimv);
  HC.topRows(dimaf) = Qafqv;
  HC.bottomRows(dimc) = Cqv;
  const Eigen::MatrixXd K_ref = - Ginv * HC;
  ASSERT_EQ(K_ref.rows(), dimaf+dimc);
  ASSERT_EQ(K_ref.cols(), 2*dimv);
  EXPECT_TRUE(K_ref.topLeftCorner(dimv, dimv).isApprox(gain.Kaq()));
  EXPECT_TRUE(K_ref.topRightCorner(dimv, dimv).isApprox(gain.Kav()));
  EXPECT_TRUE(K_ref.block(dimv, 0, dimf, dimv).isApprox(gain.Kfq()));
  EXPECT_TRUE(K_ref.block(dimv, dimv, dimf, dimv).isApprox(gain.Kfv()));
  EXPECT_TRUE(K_ref.bottomLeftCorner(dimc, dimv).isApprox(gain.Kmuq()));
  EXPECT_TRUE(K_ref.bottomRightCorner(dimc, dimv).isApprox(gain.Kmuv()));
  const Eigen::VectorXd laf = Eigen::VectorXd::Random(dimaf);
  const Eigen::VectorXd C = Eigen::VectorXd::Random(dimc);
  gain.computeFeedforward(Ginv, laf, C);
  Eigen::VectorXd h = Eigen::VectorXd::Random(dimaf+dimc);
  h.head(dimaf) = laf;
  h.tail(dimc) = C;
  const Eigen::VectorXd k_ref = - Ginv * h;
  ASSERT_EQ(k_ref.size(), dimaf+dimc);
  EXPECT_TRUE(k_ref.head(dimv).isApprox(gain.ka()));
  EXPECT_TRUE(k_ref.segment(dimv, dimf).isApprox(gain.kf()));
  EXPECT_TRUE(k_ref.tail(dimc).isApprox(gain.kmu()));
}


} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}