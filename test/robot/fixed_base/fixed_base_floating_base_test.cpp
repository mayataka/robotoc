#include <string>
#include <random>
#include <utility>
#include <vector>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "idocp/robot/floating_base.hpp"


namespace idocp {

class FixedBaseFloatingBaseTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../../urdf/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    u_ = Eigen::VectorXd::Random(dimv_);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_;
  Eigen::VectorXd u_;
};


TEST_F(FixedBaseFloatingBaseTest, constructor) {
  FloatingBase floating_base(model_);
  EXPECT_EQ(floating_base.dim_passive(), 0);
  EXPECT_TRUE(floating_base.passive_joint_indices().empty());
  EXPECT_FALSE(floating_base.has_floating_base());
}


TEST_F(FixedBaseFloatingBaseTest, moveAssign) {
  FloatingBase floating_base(model_);
  FloatingBase floating_base_ref;
  floating_base_ref = std::move(floating_base);
  EXPECT_EQ(floating_base_ref.dim_passive(), 0);
  EXPECT_TRUE(floating_base_ref.passive_joint_indices().empty());
  EXPECT_FALSE(floating_base_ref.has_floating_base());
}


TEST_F(FixedBaseFloatingBaseTest, moveConstructor) {
  FloatingBase floating_base(model_);
  FloatingBase floating_base_ref(std::move(floating_base));
  EXPECT_EQ(floating_base_ref.dim_passive(), 0);
  EXPECT_TRUE(floating_base_ref.passive_joint_indices().empty());
  EXPECT_FALSE(floating_base_ref.has_floating_base());
}


TEST_F(FixedBaseFloatingBaseTest, setPassiveTorques) {
  FloatingBase floating_base(model_);
  Eigen::VectorXd u_ref = u_;
  floating_base.setPassiveTorques(u_);
  EXPECT_TRUE(u_.isApprox(u_));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}