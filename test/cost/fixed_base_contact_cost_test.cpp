#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/contact_cost.hpp"


namespace idocp {

class FixedBaseContactCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    contact_frames_ = {18};
    robot_ = Robot(urdf_, contact_frames_, 0, 0);
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dtau_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
};


TEST_F(FixedBaseContactCostTest, zeroRefernceConstructor) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Zero(dimf);
  ContactCost cost(robot_, f_weight);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  active_contacts = {true};
  robot_.setContactStatus(active_contacts);
  ASSERT_EQ(robot_.dimf(), 3);
  const double l_ref = 0.5 * dtau_ * (f_weight.array()* (f-f_ref).array()*(f-f_ref).array()).sum();
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  lf_ref.array() = dtau_ * f_weight.array() * (f.array()-f_ref.array());
  EXPECT_TRUE(lf.isApprox(lf_ref));
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
}


TEST_F(FixedBaseContactCostTest, withRefernceConstructor) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, f_ref, f_weight);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  active_contacts = {true};
  robot_.setContactStatus(active_contacts);
  ASSERT_EQ(robot_.dimf(), 3);
  const double l_ref = 0.5 * dtau_ * (f_weight.array()* (f-f_ref).array()*(f-f_ref).array()).sum();
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  lf_ref.array() = dtau_ * f_weight.array() * (f.array()-f_ref.array());
  EXPECT_TRUE(lf.isApprox(lf_ref));
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
}


TEST_F(FixedBaseContactCostTest, setReference) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, Eigen::VectorXd::Zero(dimf), f_weight);
  cost.set_f_ref(f_ref);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  active_contacts = {true};
  robot_.setContactStatus(active_contacts);
  ASSERT_EQ(robot_.dimf(), 3);
  const double l_ref = 0.5 * dtau_ * (f_weight.array()* (f-f_ref).array()*(f-f_ref).array()).sum();
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  lf_ref.array() = dtau_ * f_weight.array() * (f.array()-f_ref.array());
  EXPECT_TRUE(lf.isApprox(lf_ref));
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
}


TEST_F(FixedBaseContactCostTest, setWeights) {
  const int dimf = robot_.max_dimf();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(dimf);
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(dimf);
  ContactCost cost(robot_, Eigen::VectorXd::Zero(dimf), Eigen::VectorXd::Zero(dimf));
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  Eigen::MatrixXd f_weight_mat = Eigen::MatrixXd::Zero(dimf, dimf);
  for (int i=0; i<dimf; ++i) {
    f_weight_mat(i, i) = f_weight(i);
  }
  const Eigen::VectorXd f = Eigen::VectorXd::Random(dimf);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(dimf);
  cost.lf(robot_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(dimf, dimf);
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  active_contacts = {true};
  robot_.setContactStatus(active_contacts);
  ASSERT_EQ(robot_.dimf(), 3);
  const double l_ref = 0.5 * dtau_ * (f_weight.array()* (f-f_ref).array()*(f-f_ref).array()).sum();
  EXPECT_FLOAT_EQ(cost.l(robot_, dtau_, f), l_ref);
  cost.lf(robot_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(dimf);
  lf_ref.array() = dtau_ * f_weight.array() * (f.array()-f_ref.array());
  EXPECT_TRUE(lf.isApprox(lf_ref));
  cost.lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
  lff.setZero();
  cost.augment_lff(robot_, dtau_, lff);
  EXPECT_TRUE(lff.isApprox(dtau_*f_weight_mat));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}