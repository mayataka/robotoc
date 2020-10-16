#include <string>
#include <random>
#include <utility>
#include <vector>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/joint_torques_lower_limit.hpp"
#include "idocp/constraints/joint_torques_upper_limit.hpp"
#include "idocp/constraints/joint_acceleration_lower_limit.hpp"
#include "idocp/constraints/joint_acceleration_upper_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

namespace idocp {

class ConstraintsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf_ = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    barrier_ = 1.0e-04;
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    amin_fixed = Eigen::VectorXd::Constant(fixed_base_robot_.dimv(), -10);
    amax_fixed = Eigen::VectorXd::Constant(fixed_base_robot_.dimv(), 10);
    amin_floating = Eigen::VectorXd::Constant(floating_base_robot_.dimv()-floating_base_robot_.dim_passive(), -10);
  }

  virtual void TearDown() {
  }

  double barrier_, dtau_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
  Eigen::VectorXd amin_fixed, amax_fixed, amin_floating;
};


TEST_F(ConstraintsTest, timestage0) {
  auto constraints = std::make_shared<Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(fixed_base_robot_);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(fixed_base_robot_);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(fixed_base_robot_);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(fixed_base_robot_);
  auto joint_torques_lower = std::make_shared<idocp::JointTorquesLowerLimit>(fixed_base_robot_);
  auto joint_torques_upper = std::make_shared<idocp::JointTorquesUpperLimit>(fixed_base_robot_);
  auto joint_accel_lower = std::make_shared<idocp::JointAccelerationLowerLimit>(fixed_base_robot_, amin_fixed);
  auto joint_accel_upper = std::make_shared<idocp::JointAccelerationUpperLimit>(fixed_base_robot_, amax_fixed);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(joint_accel_lower);
  constraints->push_back(joint_accel_upper);
  const int time_stage = 0;
  auto data = constraints->createConstraintsData(fixed_base_robot_, time_stage);
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_TRUE(data.velocity_level_data.empty());
  EXPECT_FALSE(data.acceleration_level_data.empty());
  EXPECT_EQ(data.acceleration_level_data.size(), 4);
  SplitSolution s = SplitSolution::Random(fixed_base_robot_);
  SplitDirection d = SplitDirection::Random(fixed_base_robot_);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  constraints->setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s, kkt_residual);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s.u, kkt_residual.lu);
  EXPECT_TRUE(kkt_residual.lq().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  EXPECT_TRUE(kkt_matrix.Qqq().isZero());
  EXPECT_TRUE(kkt_matrix.Qvv().isZero());
  EXPECT_FALSE(kkt_matrix.Qaa().isZero());
  EXPECT_FALSE(kkt_matrix.Quu.isZero());
  EXPECT_TRUE(kkt_residual.lq().isZero());
  EXPECT_TRUE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
}


TEST_F(ConstraintsTest, timestage1) {
  auto constraints = std::make_shared<Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(fixed_base_robot_);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(fixed_base_robot_);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(fixed_base_robot_);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(fixed_base_robot_);
  auto joint_torques_lower = std::make_shared<idocp::JointTorquesLowerLimit>(fixed_base_robot_);
  auto joint_torques_upper = std::make_shared<idocp::JointTorquesUpperLimit>(fixed_base_robot_);
  auto joint_accel_lower = std::make_shared<idocp::JointAccelerationLowerLimit>(fixed_base_robot_, amin_fixed);
  auto joint_accel_upper = std::make_shared<idocp::JointAccelerationUpperLimit>(fixed_base_robot_, amax_fixed);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(joint_accel_lower);
  constraints->push_back(joint_accel_upper);
  const int time_stage = 1;
  auto data = constraints->createConstraintsData(fixed_base_robot_, time_stage);
  EXPECT_TRUE(data.position_level_data.empty());
  EXPECT_FALSE(data.velocity_level_data.empty());
  EXPECT_FALSE(data.acceleration_level_data.empty());
  EXPECT_EQ(data.velocity_level_data.size(), 2);
  EXPECT_EQ(data.acceleration_level_data.size(), 4);
  SplitSolution s = SplitSolution::Random(fixed_base_robot_);
  SplitDirection d = SplitDirection::Random(fixed_base_robot_);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  constraints->setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s, kkt_residual);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s.u, kkt_residual.lu);
  EXPECT_TRUE(kkt_residual.lq().isZero());
  EXPECT_FALSE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  EXPECT_TRUE(kkt_matrix.Qqq().isZero());
  EXPECT_FALSE(kkt_matrix.Qvv().isZero());
  EXPECT_FALSE(kkt_matrix.Qaa().isZero());
  EXPECT_FALSE(kkt_matrix.Quu.isZero());
  EXPECT_TRUE(kkt_residual.lq().isZero());
  EXPECT_FALSE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
}


TEST_F(ConstraintsTest, timestage2) {
  auto constraints = std::make_shared<Constraints>();
  auto joint_position_lower = std::make_shared<idocp::JointPositionLowerLimit>(fixed_base_robot_);
  auto joint_position_upper = std::make_shared<idocp::JointPositionUpperLimit>(fixed_base_robot_);
  auto joint_velocity_lower = std::make_shared<idocp::JointVelocityLowerLimit>(fixed_base_robot_);
  auto joint_velocity_upper = std::make_shared<idocp::JointVelocityUpperLimit>(fixed_base_robot_);
  auto joint_torques_lower = std::make_shared<idocp::JointTorquesLowerLimit>(fixed_base_robot_);
  auto joint_torques_upper = std::make_shared<idocp::JointTorquesUpperLimit>(fixed_base_robot_);
  auto joint_accel_lower = std::make_shared<idocp::JointAccelerationLowerLimit>(fixed_base_robot_, amin_fixed);
  auto joint_accel_upper = std::make_shared<idocp::JointAccelerationUpperLimit>(fixed_base_robot_, amax_fixed);
  constraints->push_back(joint_position_lower);
  constraints->push_back(joint_position_upper);
  constraints->push_back(joint_velocity_lower);
  constraints->push_back(joint_velocity_upper);
  constraints->push_back(joint_torques_lower);
  constraints->push_back(joint_torques_upper);
  constraints->push_back(joint_accel_lower);
  constraints->push_back(joint_accel_upper);
  const int time_stage = 2;
  auto data = constraints->createConstraintsData(fixed_base_robot_, time_stage);
  EXPECT_FALSE(data.position_level_data.empty());
  EXPECT_FALSE(data.velocity_level_data.empty());
  EXPECT_FALSE(data.acceleration_level_data.empty());
  EXPECT_EQ(data.position_level_data.size(), 2);
  EXPECT_EQ(data.velocity_level_data.size(), 2);
  EXPECT_EQ(data.acceleration_level_data.size(), 4);
  SplitSolution s = SplitSolution::Random(fixed_base_robot_);
  SplitDirection d = SplitDirection::Random(fixed_base_robot_);
  KKTMatrix kkt_matrix(fixed_base_robot_);
  KKTResidual kkt_residual(fixed_base_robot_);
  constraints->setSlackAndDual(fixed_base_robot_, data, dtau_, s);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s, kkt_residual);
  constraints->augmentDualResidual(fixed_base_robot_, data, dtau_, s.u, kkt_residual.lu);
  EXPECT_FALSE(kkt_residual.lq().isZero());
  EXPECT_FALSE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s, kkt_matrix, kkt_residual);
  constraints->condenseSlackAndDual(fixed_base_robot_, data, dtau_, s.u, kkt_matrix.Quu, kkt_residual.lu);
  EXPECT_FALSE(kkt_matrix.Qqq().isZero());
  EXPECT_FALSE(kkt_matrix.Qvv().isZero());
  EXPECT_FALSE(kkt_matrix.Qaa().isZero());
  EXPECT_FALSE(kkt_matrix.Quu.isZero());
  EXPECT_FALSE(kkt_residual.lq().isZero());
  EXPECT_FALSE(kkt_residual.lv().isZero());
  EXPECT_FALSE(kkt_residual.la().isZero());
  EXPECT_FALSE(kkt_residual.lu.isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}