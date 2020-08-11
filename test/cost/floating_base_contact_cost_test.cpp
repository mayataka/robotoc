#include <string>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/contact_cost.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class FloatingBaseContactCostTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../urdf/anymal/anymal.urdf";
    contact_frames_ = {14, 24, 34, 44};
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


TEST_F(FloatingBaseContactCostTest, setWeights) {
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
  std::random_device rnd;
  std::vector<bool> active_contacts;
  for (int i=0; i<contact_frames_.size(); ++i) {
    active_contacts.push_back(rnd()%2==0);
  }
  robot_.setContactStatus(active_contacts);
  double l_ref = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        l_ref += f_weight.coeff(3*i+j) * (f.coeff(3*i+j)-f_ref.coeff(3*i+j)) 
                                       * (f.coeff(3*i+j)-f_ref.coeff(3*i+j));
      }
    }
  }
  EXPECT_DOUBLE_EQ(cost.l(robot_, data_, t_, dtau_, q, v, a, f, u), 0.5*dtau_*l_ref);
  Eigen::VectorXd lf_ref = Eigen::VectorXd::Zero(robot_.max_dimf());
  int dimf = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lf_ref.coeffRef(dimf+j) = dtau_ * f_weight.coeff(3*i+j) 
                                        * (f.coeff(3*i+j)-f_ref.coeff(3*i+j));
      }
      dimf += 3;
    }
  }
  cost.lf(robot_, data_, t_, dtau_, f, lf);
  EXPECT_TRUE(lf.isApprox(lf_ref));
  Eigen::MatrixXd lff_ref = Eigen::MatrixXd::Zero(robot_.max_dimf(), robot_.max_dimf());
  dimf = 0;
  for (int i=0; i<robot_.max_point_contacts(); ++i) {
    if (robot_.is_contact_active(i)) {
      for (int j=0; j<3; ++j) {
        lff_ref.coeffRef(dimf+j, dimf+j) = dtau_ * f_weight.coeff(3*i+j);
      }
      dimf += 3;
    }
  }
  cost.lff(robot_, data_, t_, dtau_, f, lff);
  EXPECT_TRUE(lff.isApprox(lff_ref));
  cost.augment_lff(robot_, data_, t_, dtau_, f, lff);
  EXPECT_TRUE(lff.isApprox(2*lff_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}