#include <string>
#include <random>
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
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"

#include "robotoc/robot/point_contact.hpp"

#include "urdf_factory.hpp"


namespace robotoc {

class PointContactTest : public ::testing::TestWithParam<std::pair<bool, ContactModelInfo>> {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
  }

  virtual void TearDown() {
  }

  pinocchio::Model getModel(const bool floating_base) {
    pinocchio::Model model;
    if (floating_base) {
      pinocchio::urdf::buildModel(testhelper::QuadrupedURDF(), 
                                  pinocchio::JointModelFreeFlyer(), model);
    }
    else {
      pinocchio::urdf::buildModel(testhelper::RobotManipulatorURDF(), model);
    }
    return model;
  }
};


TEST_P(PointContactTest, defaultConstructor) {
  PointContact contact;
  EXPECT_EQ(contact.contactFrameId(), 0);
  EXPECT_EQ(contact.parentJointId(), 0);
}


TEST_P(PointContactTest, constructor) {
  const bool floating_base = GetParam().first;
  const auto model = getModel(floating_base);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);
  EXPECT_EQ(contact.contactFrameId(), model.getFrameId(contact_model_info.frame));
  EXPECT_EQ(contact.parentJointId(), model.frames[model.getFrameId(contact_model_info.frame)].parent);
  EXPECT_NO_THROW(
    std::cout << contact << std::endl;
  );
}


TEST_P(PointContactTest, computeJointForcesFromContactForce) {
  const bool floating_base = GetParam().first;
  const auto model = getModel(floating_base);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(model.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  const Eigen::Vector3d fext = Eigen::Vector3d::Random();
  contact.computeJointForceFromContactForce(fext, fjoint);
  const int parent_joint_id = contact.parentJointId();
  fjoint_ref[parent_joint_id] 
      = model.frames[contact.contactFrameId()].placement.act(pinocchio::Force(fext, Eigen::Vector3d::Zero()));
  for (int i=0; i<fjoint.size(); ++i) {
    EXPECT_TRUE(fjoint[i].isApprox(fjoint_ref[i]));
  }
}


TEST_P(PointContactTest, computeBaumgarteResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero(); residual_ref.setZero();
  const Eigen::Vector3d desired_contact_position = Eigen::Vector3d::Random();
  contact.computeBaumgarteResidual(model, data, desired_contact_position, residual);
  residual_ref 
      = pinocchio::getFrameClassicalAcceleration(model, data, contact.contactFrameId(), pinocchio::LOCAL).linear()
          + contact_model_info.baumgarte_velocity_gain 
              * pinocchio::getFrameVelocity(model, data, contact.contactFrameId(), pinocchio::LOCAL).linear()
          + contact_model_info.baumgarte_position_gain
              * (data.oMf[contact.contactFrameId()].translation()-desired_contact_position);
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(PointContactTest, computeBaumgarteDerivatives) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeBaumgarteDerivatives(model, data, baum_partial_dq, baum_partial_dv, baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_da = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameAccelerationDerivatives(model, data, contact.contactFrameId(), pinocchio::LOCAL,
                                             frame_v_partial_dq, frame_a_partial_dq, 
                                             frame_a_partial_dv, frame_a_partial_da);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact.contactFrameId(), pinocchio::LOCAL, J_frame);
  pinocchio::Motion v_frame = pinocchio::getFrameVelocity(model, data, contact.contactFrameId(), pinocchio::LOCAL);
  Eigen::Matrix3d v_linear_skew, v_angular_skew;
  v_linear_skew.setZero(); v_angular_skew.setZero();
  pinocchio::skew(v_frame.linear(), v_linear_skew);
  pinocchio::skew(v_frame.angular(), v_angular_skew);
  baum_partial_dq_ref 
      = frame_a_partial_dq.template topRows<3>()
          + v_angular_skew * frame_v_partial_dq.template topRows<3>()
          - v_linear_skew * frame_v_partial_dq.template bottomRows<3>()
          + contact_model_info.baumgarte_velocity_gain * frame_v_partial_dq.template topRows<3>()
          + contact_model_info.baumgarte_position_gain * data.oMf[contact.contactFrameId()].rotation() * J_frame.template topRows<3>();
  baum_partial_dv_ref 
      = frame_a_partial_dv.template topRows<3>()
          + v_angular_skew * J_frame.template topRows<3>()
          - v_linear_skew * J_frame.template bottomRows<3>()
          + contact_model_info.baumgarte_velocity_gain * frame_a_partial_da.template topRows<3>();
  baum_partial_da_ref 
      = frame_a_partial_da.template topRows<3>();
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


TEST_P(PointContactTest, computeContactVelocityResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero(); residual_ref.setZero();
  contact.computeContactVelocityResidual(model, data, residual);
  residual_ref = pinocchio::getFrameVelocity(model, data, contact.contactFrameId(), pinocchio::LOCAL).linear();
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(PointContactTest, computeContactVelocityDerivatives) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  Eigen::MatrixXd vel_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd vel_partial_dv = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeContactVelocityDerivatives(model, data, vel_partial_dq, vel_partial_dv);
  Eigen::MatrixXd vel_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd vel_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameVelocityDerivatives(model, data, contact.contactFrameId(), pinocchio::LOCAL,
                                         frame_v_partial_dq, frame_v_partial_dv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  vel_partial_dq_ref = frame_v_partial_dq.template topRows<3>();
  vel_partial_dv_ref = frame_v_partial_dv.template topRows<3>();
  EXPECT_TRUE(vel_partial_dq_ref.isApprox(vel_partial_dq));
  EXPECT_TRUE(vel_partial_dv_ref.isApprox(vel_partial_dv));
}


TEST_P(PointContactTest, computeContactPositionResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const Eigen::Vector3d desired_contact_position = Eigen::Vector3d::Random();
  contact.computeContactPositionResidual(model, data, desired_contact_position, residual);
  residual_ref = (data.oMf[contact.contactFrameId()].translation()-desired_contact_position);
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(PointContactTest, computeContactPositionDerivative) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  PointContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  Eigen::MatrixXd position_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeContactPositionDerivative(model, data, position_partial_dq);
  Eigen::MatrixXd position_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact.contactFrameId(), pinocchio::LOCAL, J_frame);
  position_partial_dq_ref = data.oMf[contact.contactFrameId()].rotation() * J_frame.template topRows<3>();
  EXPECT_TRUE(position_partial_dq_ref.isApprox(position_partial_dq));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, PointContactTest, 
  ::testing::Values(std::make_pair(false, ContactModelInfo("iiwa_link_ee_kuka", 0.001)),
                    std::make_pair(true, ContactModelInfo("LF_FOOT", 0.001)))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}