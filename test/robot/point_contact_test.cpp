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

#include "idocp/robot/point_contact.hpp"


namespace idocp {

class PointContactTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(fixed_base_urdf, fixed_base_robot);
    pinocchio::urdf::buildModel(floating_base_urdf, floating_base_robot);
    fixed_base_data = pinocchio::Data(fixed_base_robot);
    floating_base_data = pinocchio::Data(floating_base_robot);
    fixed_base_contact_frames = {18};
    floating_base_contact_frames = {14, 24, 34, 44};
    baumgarte_weight_on_velocity = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
    baumgarte_weight_on_position = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testConstructor(pinocchio::Model& model, pinocchio::Data& data, 
                       const int contact_frame_id) const;
  void testComputeJointForcesFromContactForce(pinocchio::Model& model, 
                                              pinocchio::Data& data, 
                                              const int contact_frame_id) const;
  void testGetContactJacobian(pinocchio::Model& model, pinocchio::Data& data, 
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


TEST_F(PointContactTest, defaultConstructor) {
  PointContact contact;
  EXPECT_EQ(contact.contact_frame_id(), 0);
  EXPECT_EQ(contact.parent_joint_id(), 0);
}


void PointContactTest::testConstructor(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id);
  EXPECT_EQ(contact.parent_joint_id(), model.frames[contact_frame_id].parent);
}


void PointContactTest::testComputeJointForcesFromContactForce(
    pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(
            model.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  const Eigen::Vector3d fext = Eigen::Vector3d::Random();
  contact.computeJointForceFromContactForce(fext, fjoint);
  const int parent_joint_id = contact.parent_joint_id();
  fjoint_ref[parent_joint_id] 
      = model.frames[contact_frame_id].placement.act(
          pinocchio::Force(fext, Eigen::Vector3d::Zero()));
  for (int i=0; i<fjoint.size(); ++i) {
    EXPECT_TRUE(fjoint[i].isApprox(fjoint_ref[i]));
  }
}


void PointContactTest::testGetContactJacobian(
    pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q);
  pinocchio::computeJointJacobians(model, data, q);
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd J_trans = Eigen::MatrixXd::Zero(dimv, 3);
  Eigen::MatrixXd J_ref = Eigen::MatrixXd::Zero(6, dimv);
  contact.getContactJacobian(model, data, J);
  pinocchio::getFrameJacobian(model, data, contact_frame_id, pinocchio::LOCAL, J_ref);
  EXPECT_TRUE(J.isApprox(J_ref.topRows(3)));
  const bool transpose = true;
  contact.getContactJacobian(model, data, J_trans, transpose);
  EXPECT_TRUE(J_trans.isApprox(J_ref.topRows(3).transpose()));
  const int block_rows_begin = rand() % 5;
  const int block_cols_begin = rand() % 5;
  Eigen::MatrixXd J_block 
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv+2*block_cols_begin);
  contact.getContactJacobian(model, data, J_block.block(block_rows_begin, 
                                                         block_cols_begin,
                                                         3, dimv));
  EXPECT_TRUE(J_block.block(block_rows_begin, block_cols_begin, 3, dimv).
              isApprox(J));
  Eigen::MatrixXd J_trans_block 
      = Eigen::MatrixXd::Zero(dimv+2*block_cols_begin, 3+2*block_rows_begin);
  contact.getContactJacobian(model, data, J_trans_block.block(block_cols_begin, 
                                                               block_rows_begin,
                                                               dimv, 3), transpose);
  EXPECT_TRUE(J_block.isApprox(J_trans_block.transpose()));
  J_block.block(block_rows_begin, block_cols_begin, 3, dimv) -= J;
  EXPECT_TRUE(J_block.isZero());
}


void PointContactTest::testBaumgarteResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const Eigen::Vector3d contact_point = Eigen::Vector3d::Random();
  contact.computeBaumgarteResidual(model, data, contact_point, residual);
  residual_ref 
      = pinocchio::getFrameClassicalAcceleration(model, data, contact_frame_id, 
                                                 pinocchio::LOCAL).linear()
          + baumgarte_weight_on_velocity
              * pinocchio::getFrameVelocity(model, data, contact_frame_id, 
                                            pinocchio::LOCAL).linear()
          + baumgarte_weight_on_position 
              * (data.oMf[contact_frame_id].translation()
                 -contact_point);
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeBaumgarteResidual(model, data, contact_point, residuals.segment<3>(5));
  EXPECT_TRUE(residuals.head(5).isZero());
  EXPECT_TRUE(residuals.segment<3>(5).isApprox(residual_ref));
}


void PointContactTest::testBaumgarteDerivative(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeBaumgarteDerivatives(model, data, baum_partial_dq, 
                                      baum_partial_dv, baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(3, dimv);
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
  pinocchio::Motion v_frame = pinocchio::getFrameVelocity(model, data, 
                                                          contact_frame_id, 
                                                          pinocchio::LOCAL);
  Eigen::Matrix3d v_linear_skew, v_angular_skew;
  v_linear_skew.setZero();
  v_angular_skew.setZero();
  pinocchio::skew(v_frame.linear(), v_linear_skew);
  pinocchio::skew(v_frame.angular(), v_angular_skew);
  baum_partial_dq_ref 
      = frame_a_partial_dq.template topRows<3>()
          + v_angular_skew * frame_v_partial_dq.template topRows<3>()
          - v_linear_skew * frame_v_partial_dq.template bottomRows<3>()
          + baumgarte_weight_on_velocity 
              * frame_v_partial_dq.template topRows<3>()
          + baumgarte_weight_on_position 
              * data.oMf[contact_frame_id].rotation()
              * J_frame.template topRows<3>();
  baum_partial_dv_ref 
      = frame_a_partial_dv.template topRows<3>()
          + v_angular_skew * J_frame.template topRows<3>()
          - v_linear_skew * J_frame.template bottomRows<3>()
          + baumgarte_weight_on_velocity
              * frame_a_partial_da.template topRows<3>();
  baum_partial_da_ref 
      = frame_a_partial_da.template topRows<3>();
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


void PointContactTest::testContactVelocityResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  contact.computeContactVelocityResidual(model, data, residual);
  residual_ref = pinocchio::getFrameVelocity(model, data, contact_frame_id, 
                                             pinocchio::LOCAL).linear();
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeContactVelocityResidual(model, data, residuals.segment<3>(5));
  EXPECT_TRUE(residuals.head(5).isZero());
  EXPECT_TRUE(residuals.segment<3>(5).isApprox(residual_ref));
  EXPECT_TRUE(residuals.tail(2).isZero());
}


void PointContactTest::testContactVelocityDerivatives(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd vel_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd vel_partial_dv = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeContactVelocityDerivatives(model, data, 
                                            vel_partial_dq, vel_partial_dv);
  Eigen::MatrixXd vel_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd vel_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv);
  Eigen::MatrixXd frame_v_partial_dv = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameVelocityDerivatives(model, data, contact_frame_id, 
                                         pinocchio::LOCAL,
                                         frame_v_partial_dq, 
                                         frame_v_partial_dv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  vel_partial_dq_ref = frame_v_partial_dq.template topRows<3>();
  vel_partial_dv_ref = frame_v_partial_dv.template topRows<3>();
  EXPECT_TRUE(vel_partial_dq_ref.isApprox(vel_partial_dq));
  EXPECT_TRUE(vel_partial_dv_ref.isApprox(vel_partial_dv));
}


void PointContactTest::testContactResidual(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  const Eigen::Vector3d contact_point = Eigen::Vector3d::Random();
  contact.computeContactResidual(model, data, contact_point, residual);
  residual_ref = (data.oMf[contact_frame_id].translation()-contact_point);
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeContactResidual(model, data, contact_point, residuals.segment<3>(5));
  EXPECT_TRUE(residuals.head(5).isZero());
  EXPECT_TRUE(residuals.segment<3>(5).isApprox(residual_ref));
  EXPECT_TRUE(residuals.tail(2).isZero());
}


void PointContactTest::testContactDerivatives(pinocchio::Model& model, pinocchio::Data& data, const int contact_frame_id) const {
  PointContact contact(model, contact_frame_id, baumgarte_weight_on_velocity, baumgarte_weight_on_position);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
  const int dimv = model.nv;
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  pinocchio::updateFramePlacement(model, data, contact_frame_id);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd contact_partial_dq = Eigen::MatrixXd::Zero(3, dimv);
  contact.computeContactDerivative(model, data, contact_partial_dq);
  Eigen::MatrixXd contact_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv);
  pinocchio::getFrameJacobian(model, data, contact_frame_id, 
                              pinocchio::LOCAL, J_frame);
  contact_partial_dq_ref 
      = data.oMf[contact_frame_id].rotation() * J_frame.template topRows<3>();
  EXPECT_TRUE(contact_partial_dq_ref.isApprox(contact_partial_dq));
}


TEST_F(PointContactTest, test) {
  for (const auto frame : fixed_base_contact_frames) {
    pinocchio::Model robot = fixed_base_robot;
    pinocchio::Data data = fixed_base_data;
    testConstructor(robot, data, frame);
    testComputeJointForcesFromContactForce(robot, data, frame);
    testGetContactJacobian(robot, data, frame);
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
    testComputeJointForcesFromContactForce(robot, data, frame);
    testGetContactJacobian(robot, data, frame);
    testBaumgarteResidual(robot, data, frame);
    testBaumgarteDerivative(robot, data, frame);
    testContactVelocityResidual(robot, data, frame);
    testContactVelocityDerivatives(robot, data, frame);
    testContactResidual(robot, data, frame);
    testContactDerivatives(robot, data, frame);
  }
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}