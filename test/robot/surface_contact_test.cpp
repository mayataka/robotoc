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

#include "robotoc/robot/surface_contact.hpp"

#include "urdf_factory.hpp"


namespace robotoc {

class SurfaceContactTest : public ::testing::TestWithParam<std::pair<bool, ContactModelInfo>> {
protected:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

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


TEST_P(SurfaceContactTest, defaultConstructor) {
  SurfaceContact contact;
  EXPECT_EQ(contact.contactFrameId(), 0);
  EXPECT_EQ(contact.parentJointId(), 0);
}


TEST_P(SurfaceContactTest, constructor) {
  const bool floating_base = GetParam().first;
  const auto model = getModel(floating_base);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  EXPECT_EQ(contact.contactFrameId(), model.getFrameId(contact_model_info.frame));
  EXPECT_EQ(contact.parentJointId(), model.frames[model.getFrameId(contact_model_info.frame)].parent);
  EXPECT_NO_THROW(
    std::cout << contact << std::endl;
  );
}


TEST_P(SurfaceContactTest, computeJointForcesFromContactWrench) {
  const bool floating_base = GetParam().first;
  const auto model = getModel(floating_base);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(model.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  const Vector6d fext = Vector6d::Random();
  contact.computeJointForceFromContactWrench(fext, fjoint);
  const int parent_joint_id = contact.parentJointId();
  fjoint_ref[parent_joint_id] 
      = model.frames[contact.contactFrameId()].placement.act(pinocchio::Force(fext));
  for (int i=0; i<fjoint.size(); ++i) {
    EXPECT_TRUE(fjoint[i].isApprox(fjoint_ref[i]));
  }
}


TEST_P(SurfaceContactTest, computeBaumgarteResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Vector6d residual, residual_ref;
  residual.setZero(); residual_ref.setZero();
  const SE3 desired_contact_placement = SE3::Random();
  contact.computeBaumgarteResidual(model, data, desired_contact_placement, residual);
  residual_ref = pinocchio::getFrameAcceleration(model, data, contact.contactFrameId(), pinocchio::LOCAL).toVector()
                  + contact_model_info.baumgarte_velocity_gain
                      * pinocchio::getFrameVelocity(model, data, contact.contactFrameId(), pinocchio::LOCAL).toVector()
                  + contact_model_info.baumgarte_position_gain 
                      * pinocchio::log6(desired_contact_placement.inverse()*data.oMf[contact.contactFrameId()]).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(SurfaceContactTest, computeBaumgarteDerivatives) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(6, dimv);
  const SE3 desired_contact_placement = SE3::Random();
  Vector6d residual;
  contact.computeBaumgarteResidual(model, data, desired_contact_placement, residual);
  contact.computeBaumgarteDerivatives(model, data, baum_partial_dq, baum_partial_dv, baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_da = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameAccelerationDerivatives(model, data, contact.contactFrameId(), pinocchio::LOCAL,
                                             frame_v_partial_dq, frame_a_partial_dq, 
                                             frame_a_partial_dv, frame_a_partial_da);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact.contactFrameId(), pinocchio::LOCAL, J_frame);
  Eigen::MatrixXd Jlog6(Eigen::MatrixXd::Zero(6, 6));
  pinocchio::Jlog6(desired_contact_placement.inverse()*data.oMf[contact.contactFrameId()], Jlog6);
  baum_partial_dq_ref 
      = frame_a_partial_dq
          + contact_model_info.baumgarte_velocity_gain * frame_v_partial_dq
          + contact_model_info.baumgarte_position_gain * Jlog6 * J_frame;
  baum_partial_dv_ref 
      = frame_a_partial_dv
          + contact_model_info.baumgarte_velocity_gain * frame_a_partial_da;
  baum_partial_da_ref = frame_a_partial_da;
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


TEST_P(SurfaceContactTest, computeContactVelocityResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Vector6d residual, residual_ref;
  residual.setZero(); residual_ref.setZero();
  contact.computeContactVelocityResidual(model, data, residual);
  residual_ref = pinocchio::getFrameVelocity(model, data, contact.contactFrameId(), pinocchio::LOCAL).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(SurfaceContactTest, computeContactVelocityDerivatives) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  Eigen::MatrixXd vel_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd vel_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  contact.computeContactVelocityDerivatives(model, data, vel_partial_dq, vel_partial_dv);
  Eigen::MatrixXd vel_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd vel_partial_dv_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameVelocityDerivatives(model, data, contact.contactFrameId(), pinocchio::LOCAL,
                                         frame_v_partial_dq, frame_v_partial_dv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  vel_partial_dq_ref = frame_v_partial_dq;
  vel_partial_dv_ref = frame_v_partial_dv;
  EXPECT_TRUE(vel_partial_dq_ref.isApprox(vel_partial_dq));
  EXPECT_TRUE(vel_partial_dv_ref.isApprox(vel_partial_dv));
}


TEST_P(SurfaceContactTest, computeContactPositionResidual) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parentJointId();
  Vector6d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const SE3 desired_contact_placement = SE3::Random();
  contact.computeContactPositionResidual(model, data, desired_contact_placement, residual);
  residual_ref = pinocchio::log6(desired_contact_placement.inverse()*data.oMf[contact.contactFrameId()]).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
}


TEST_P(SurfaceContactTest, computeContactPositionDerivative) {
  const bool floating_base = GetParam().first;
  auto model = getModel(floating_base);
  auto data = pinocchio::Data(model);
  const auto contact_model_info = GetParam().second;
  SurfaceContact contact(model, contact_model_info);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact.contactFrameId());
  const int parent_joint_id = contact.parentJointId();
  const SE3 desired_contact_placement = SE3::Random();
  Vector6d residual;
  contact.computeContactPositionResidual(model, data, desired_contact_placement, residual);
  Eigen::MatrixXd position_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  contact.computeContactPositionDerivative(model, data, position_partial_dq);
  Eigen::MatrixXd position_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact.contactFrameId(), pinocchio::LOCAL, J_frame);
  Eigen::MatrixXd Jlog6(Eigen::MatrixXd::Zero(6, 6));
  pinocchio::Jlog6(desired_contact_placement.inverse()*data.oMf[contact.contactFrameId()], Jlog6);
  position_partial_dq_ref = Jlog6 * J_frame;
  EXPECT_TRUE(position_partial_dq_ref.isApprox(position_partial_dq));
}


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, SurfaceContactTest, 
  ::testing::Values(std::make_pair(false, ContactModelInfo("iiwa_link_ee_kuka", 0.001)),
                    std::make_pair(true, ContactModelInfo("LF_FOOT", 0.001)))
);

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}