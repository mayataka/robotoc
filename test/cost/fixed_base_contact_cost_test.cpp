#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"


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
    t_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    data_ = CostFunctionData(robot_);
  }

  virtual void TearDown() {
  }

  double dtau_, t_;
  std::vector<int> contact_frames_;
  std::string urdf_;
  Robot robot_;
  CostFunctionData data_;
};


TEST_F(FixedBaseContactCostTest, setWeights) {
  const int dimq = robot_.dimq();
  const int dimv = robot_.dimv();
  const Eigen::VectorXd f_weight = Eigen::VectorXd::Random(robot_.max_dimf());
  const Eigen::VectorXd f_ref = Eigen::VectorXd::Random(robot_.max_dimf());
  ContactCost cost(robot_);
  cost.set_f_weight(f_weight);
  cost.set_f_ref(f_ref);
  const Eigen::VectorXd q = Eigen::VectorXd::Random(dimq);
  const Eigen::VectorXd v = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(dimv);
  const Eigen::VectorXd f = Eigen::VectorXd::Random(robot_.max_dimf());
  const Eigen::VectorXd u = Eigen::VectorXd::Random(dimv);
  ASSERT_EQ(robot_.dimf(), 0);
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, q, v, a, f, u), 0);
  Eigen::VectorXd lf = Eigen::VectorXd::Zero(robot_.max_dimf());
  cost.lf(robot_, data_, t_, dtau_, f, lf);
  EXPECT_TRUE(lf.isZero());
  Eigen::MatrixXd lff = Eigen::MatrixXd::Zero(robot_.max_dimf(), robot_.max_dimf());
  cost.lff(robot_, data_, t_, dtau_, f, lff);
  EXPECT_TRUE(lff.isZero());
  std::vector<bool> active_contacts;
  active_contacts = {true};
  robot_.setContactStatus(active_contacts);
  ASSERT_EQ(robot_.dimf(), 3);
  const double l_ref = 0.5 * dtau_ * (f_weight.array()* (f-f_ref).array()*(f-f_ref).array()).sum();
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, q, v, a, f, u), l_ref);
  cost.lf(robot_, data_, t_, dtau_, f, lf);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(robot_.max_dimf());
  lf_ref.array() = dtau_ * f_weight.asDiagonal() * (f-f_ref);
  EXPECT_TRUE(lf.isApprox(lf_ref));
  cost.lff(robot_, data_, t_, dtau_, f, lff);
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(robot_.max_dimf(), robot_.max_dimf());
  lff_ref = dtau_*f_weight.asDiagonal();
  EXPECT_TRUE(lff.isApprox(lff_ref));
  cost.augment_lff(robot_, data_, t_, dtau_, f, lff);
  EXPECT_TRUE(lff.isApprox(2*lff_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}