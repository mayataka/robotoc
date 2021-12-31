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


namespace robotoc {

class SurfaceContactTest : public ::testing::Test {
protected:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(fixed_base_urdf, fixed_base_robot);
    pinocchio::urdf::buildModel(floating_base_urdf, 
                                pinocchio::JointModelFreeFlyer(), 
                                floating_base_robot);
    fixed_base_data = pinocchio::Data(fixed_base_robot);
    floating_base_data = pinocchio::Data(floating_base_robot);
    fixed_base_contact_frames = {18};
    // floating_base_contact_frames = {12, 22, 32, 42};
    floating_base_contact_frames = {12, 22};
    baumgarte_weight_on_velocity = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
    baumgarte_weight_on_position = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testConstructor(pinocchio::Model& model, pinocchio::Data& data, 
                       const int contact_frame_id) const;
  void testComputeJointForcesFromContactWrench(pinocchio::Model& model, 
                                              pinocchio::Data& data, 
                                              const int contact_frame_id) const;
  void testBaumgarteResidual(pinocchio::Model& model, pinocchio::Data& data, 
                              const int contact_frame_id) const;
  void testBaumgarteDerivative(pinocchio::Model& model, pinocchio::Data& data, 
                               const int contact_frame_id) const;
  void testContactVelocityResidual(pinocchio::Model& model, 
                                   pinocchio::Data& data, 
                                   const int contact_frame_id) const;
  void testContactVelocityDerivatives(pinocchio::Model& model, 
                                      pinocchio::Data& data, 
                                      const int contact_frame_id) const;
  void testContactResidual(pinocchio::Model& model, pinocchio::Data& data, 
                           const int contact_frame_id) const;
  void testContactDerivatives(pinocchio::Model& model, pinocchio::Data& data, 
                              const int contact_frame_id) const;

  std::string fixed_base_urdf, floating_base_urdf;
  pinocchio::Model fixed_base_robot, floating_base_robot;
  pinocchio::Data fixed_base_data, floating_base_data;
  std::vector<int> fixed_base_contact_frames, floating_base_contact_frames;
  double baumgarte_weight_on_velocity, baumgarte_weight_on_position;
};


TEST_F(SurfaceContactTest, defaultConstructor) {
  SurfaceContact contact;
  EXPECT_EQ(contact.contact_frame_id(), 0);
  EXPECT_EQ(contact.parent_joint_id(), 0);
}


void SurfaceContactTest::testConstructor(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id);
  EXPECT_EQ(contact.parent_joint_id(), model.frames[contact_frame_id].parent);
  EXPECT_NO_THROW(
    std::cout << contact << std::endl;
  );
}


void SurfaceContactTest::testComputeJointForcesFromContactWrench(
    pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(
            model.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  const Vector6d fext = Vector6d::Random();
  contact.computeJointForceFromContactWrench(fext, fjoint);
  const int parent_joint_id = contact.parent_joint_id();
  fjoint_ref[parent_joint_id] 
      = model.frames[contact_frame_id].placement.act(pinocchio::Force(fext));
  for (int i=0; i<fjoint.size(); ++i) {
    EXPECT_TRUE(fjoint[i].isApprox(fjoint_ref[i]));
  }
}


void SurfaceContactTest::testBaumgarteResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Vector6d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const SE3 contact_placement = SE3::Random();
  contact.computeBaumgarteResidual(model, data, contact_placement, residual);
  residual_ref = pinocchio::getFrameAcceleration(model, data, contact_frame_id, 
                                                 pinocchio::LOCAL).toVector()
                  + baumgarte_weight_on_velocity
                      * pinocchio::getFrameVelocity(model, data, contact_frame_id, 
                                                    pinocchio::LOCAL).toVector()
                  + baumgarte_weight_on_position 
                      * pinocchio::log6(contact_placement.inverse()*data.oMf[contact_frame_id]).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeBaumgarteResidual(model, data, contact_placement, residuals.segment<6>(3));
  EXPECT_TRUE(residuals.head(3).isZero());
  EXPECT_TRUE(residuals.segment<6>(3).isApprox(residual_ref));
}


void SurfaceContactTest::testBaumgarteDerivative(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(6, dimv);
  const SE3 contact_placement = SE3::Random();
  Vector6d residual;
  contact.computeBaumgarteResidual(model, data, contact_placement, residual);
  contact.computeBaumgarteDerivatives(model, data, baum_partial_dq, 
                                      baum_partial_dv, baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_a_partial_da = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameAccelerationDerivatives(model, data, contact_frame_id, 
                                             pinocchio::LOCAL,
                                             frame_v_partial_dq, 
                                             frame_a_partial_dq, 
                                             frame_a_partial_dv, 
                                             frame_a_partial_da);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact_frame_id, 
                              pinocchio::LOCAL, J_frame);
  Eigen::MatrixXd Jlog6(Eigen::MatrixXd::Zero(6, 6));
  pinocchio::Jlog6(contact_placement.inverse()*data.oMf[contact_frame_id], Jlog6);
  baum_partial_dq_ref 
      = frame_a_partial_dq
          + baumgarte_weight_on_velocity * frame_v_partial_dq
          + baumgarte_weight_on_position * Jlog6 * J_frame;
  baum_partial_dv_ref 
      = frame_a_partial_dv
          + baumgarte_weight_on_velocity * frame_a_partial_da;
  baum_partial_da_ref = frame_a_partial_da;
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


void SurfaceContactTest::testContactVelocityResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Vector6d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  contact.computeContactVelocityResidual(model, data, residual);
  residual_ref = pinocchio::getFrameVelocity(model, data, contact_frame_id, 
                                             pinocchio::LOCAL).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeContactVelocityResidual(model, data, residuals.segment<6>(3));
  EXPECT_TRUE(residuals.head(3).isZero());
  EXPECT_TRUE(residuals.segment<6>(3).isApprox(residual_ref));
}


void SurfaceContactTest::testContactVelocityDerivatives(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd vel_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd vel_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  contact.computeContactVelocityDerivatives(model, data, 
                                            vel_partial_dq, vel_partial_dv);
  Eigen::MatrixXd vel_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd vel_partial_dv_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameVelocityDerivatives(model, data, contact_frame_id, 
                                         pinocchio::LOCAL,
                                         frame_v_partial_dq, 
                                         frame_v_partial_dv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  vel_partial_dq_ref = frame_v_partial_dq;
  vel_partial_dv_ref = frame_v_partial_dv;
  EXPECT_TRUE(vel_partial_dq_ref.isApprox(vel_partial_dq));
  EXPECT_TRUE(vel_partial_dv_ref.isApprox(vel_partial_dv));
}


void SurfaceContactTest::testContactResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Vector6d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const SE3 contact_placement = SE3::Random();
  contact.computeContactPositionResidual(model, data, contact_placement, residual);
  residual_ref = pinocchio::log6(contact_placement.inverse()*data.oMf[contact_frame_id]).toVector();
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeContactPositionResidual(model, data, contact_placement, residuals.segment<6>(3));
  EXPECT_TRUE(residuals.head(3).isZero());
  EXPECT_TRUE(residuals.segment<6>(3).isApprox(residual_ref));
}


void SurfaceContactTest::testContactDerivatives(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  SurfaceContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  const SE3 contact_placement = SE3::Random();
  Vector6d residual;
  contact.computeContactPositionResidual(model, data, contact_placement, residual);
  Eigen::MatrixXd position_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  contact.computeContactPositionDerivative(model, data, position_partial_dq);
  Eigen::MatrixXd position_partial_dq_ref = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact_frame_id, 
                              pinocchio::LOCAL, J_frame);
  Eigen::MatrixXd Jlog6(Eigen::MatrixXd::Zero(6, 6));
  pinocchio::Jlog6(contact_placement.inverse()*data.oMf[contact_frame_id], Jlog6);
  position_partial_dq_ref = Jlog6 * J_frame;
  EXPECT_TRUE(position_partial_dq_ref.isApprox(position_partial_dq));
}


TEST_F(SurfaceContactTest, test) {
  for (const auto frame : fixed_base_contact_frames) {
    pinocchio::Model robot = fixed_base_robot;
    pinocchio::Data data = fixed_base_data;
    testConstructor(robot, data, frame);
    testComputeJointForcesFromContactWrench(robot, data, frame);
    testBaumgarteResidual(robot, data, frame);
    testBaumgarteDerivative(robot, data, frame);
    testContactVelocityResidual(robot, data, frame);
    testContactVelocityDerivatives(robot, data, frame);
    testContactResidual(robot, data, frame);
    testContactDerivatives(robot, data, frame);
  }
  for (const auto frame : floating_base_contact_frames) {
    pinocchio::Model robot = floating_base_robot;
    pinocchio::Data data = floating_base_data;
    testConstructor(robot, data, frame);
    testComputeJointForcesFromContactWrench(robot, data, frame);
    testBaumgarteResidual(robot, data, frame);
    testBaumgarteDerivative(robot, data, frame);
    testContactVelocityResidual(robot, data, frame);
    testContactVelocityDerivatives(robot, data, frame);
    testContactResidual(robot, data, frame);
    testContactDerivatives(robot, data, frame);
  }
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}