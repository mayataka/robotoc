#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "robot/passive_joints.hpp"


namespace invdynocp {

class PassiveJointsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_file_name_ = "../../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    njoints_ = model_.njoints;
    u_ = Eigen::VectorXd::Random(dimv_);
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[1];
  }

  virtual void TearDown() {
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, njoints_, contact_frame_id_;
  double baumgarte_alpha_, baumgarte_beta_;
  Eigen::VectorXd u_;
};


TEST_F(PassiveJointsTest, passive_torque_indices) {
  PassiveJoints passive_joints_(model_);
  EXPECT_EQ(passive_joints_.dim_passive(), 6);
  for (int i=0; i<passive_joints_.dim_passive(); ++i) {
    EXPECT_EQ(passive_joints_.passive_torque_indices()[i], i);
  }
}


TEST_F(PassiveJointsTest, setPassiveTorques) {
  PassiveJoints passive_joints_(model_);
  passive_joints_.setPassiveTorques(u_);
  for (int i=0; i<passive_joints_.dim_passive(); ++i) {
    EXPECT_DOUBLE_EQ(u_[i], 0);
  }
}


TEST_F(PassiveJointsTest, computePassiveConstraintViolation) {
  PassiveJoints passive_joints_(model_);
  Eigen::VectorXd violation = Eigen::VectorXd::Zero(passive_joints_.dim_passive());
  passive_joints_.computePassiveConstraintViolation(u_, violation);
  for (int i=0; i<passive_joints_.dim_passive(); ++i) {
    EXPECT_DOUBLE_EQ(violation[i], u_[i]);
  }
}


} // namespace invdynocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}