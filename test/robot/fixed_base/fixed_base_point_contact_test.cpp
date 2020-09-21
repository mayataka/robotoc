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
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/frames-derivatives.hpp"

#include "idocp/robot/point_contact.hpp"


namespace idocp {

class FixedBasePointContactTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_ = "../../urdf/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    njoints_ = model_.njoints;
    contact_frame_id_ = 2 + (int)((rnd()%(model_.nframes-2)) / 2) * 2;
    q_ = pinocchio::randomConfiguration(model_, 
                                        -Eigen::VectorXd::Ones(dimq_), 
                                        Eigen::VectorXd::Ones(dimq_));
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    friction_coeff_ = std::abs(Eigen::VectorXd::Random(2)[1]);
    restitution_coeff_ = std::abs(Eigen::VectorXd::Random(2)[1]);
    while (restitution_coeff_ > 1) {
      restitution_coeff_ = std::abs(Eigen::VectorXd::Random(2)[1]);
    }
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, njoints_, contact_frame_id_;
  double friction_coeff_, restitution_coeff_;
  Eigen::VectorXd q_, v_, a_;
};


TEST_F(FixedBasePointContactTest, defaultConstructor) {
  PointContact contact;
  // Check the default constructor works appropriately.
  EXPECT_FALSE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), 0);
  EXPECT_EQ(contact.parent_joint_id(), 0);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), 0);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), 0);
  EXPECT_TRUE(contact.contactPoint().isApprox(Eigen::Vector3d::Zero()));
}


TEST_F(FixedBasePointContactTest, constructor) {
  PointContact contact(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  // Check the constructor works appropriately.
  EXPECT_FALSE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), friction_coeff_);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), restitution_coeff_);
  EXPECT_TRUE(contact.contactPoint().isApprox(Eigen::Vector3d::Zero()));
  // Check set parameters set appropriately.
  const double friction_coeff_tmp = std::abs(Eigen::VectorXd::Random(2)[1]);
  const double restitution_coeff_tmp = std::abs(Eigen::VectorXd::Random(2)[1]);
  const Eigen::Vector3d contact_point_tmp = Eigen::Vector3d::Random();
  contact.setFrictionCoefficient(friction_coeff_tmp);
  contact.setRestitutionCoefficient(restitution_coeff_tmp);
  contact.setContactPoint(contact_point_tmp);
  EXPECT_DOUBLE_EQ(friction_coeff_tmp, contact.frictionCoefficient());
  EXPECT_DOUBLE_EQ(restitution_coeff_tmp, contact.restitutionCoefficient());
  EXPECT_TRUE(contact_point_tmp.isApprox(contact.contactPoint()));
  // Check activate works appropriately.
  contact.activate();
  EXPECT_TRUE(contact.isActive());
  contact.deactivate();
  EXPECT_FALSE(contact.isActive());
  // Check the contact point assignment by kinematics.
  pinocchio::forwardKinematics(model_, data_, q_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  contact.setContactPointByCurrentKinematics(data_);
  EXPECT_TRUE(
      contact.contactPoint().isApprox(
          data_.oMf[contact_frame_id_].translation()));
}


TEST_F(FixedBasePointContactTest, copyConstructor) {
  PointContact contact_ref(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  contact_ref.activate();
  const Eigen::Vector3d contact_point_ref = Eigen::Vector3d::Random();
  contact_ref.setContactPoint(contact_point_ref);
  PointContact contact(contact_ref);
  // Check the constructor works well.
  EXPECT_TRUE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), friction_coeff_);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), restitution_coeff_);
  EXPECT_TRUE(contact.contactPoint().isApprox(contact_point_ref));
}


TEST_F(FixedBasePointContactTest, assign) {
  PointContact contact_ref(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  contact_ref.activate();
  const Eigen::Vector3d contact_point_ref = Eigen::Vector3d::Random();
  contact_ref.setContactPoint(contact_point_ref);
  PointContact contact;
  contact = contact_ref;
  // Check the constructor works well.
  EXPECT_TRUE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), friction_coeff_);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), restitution_coeff_);
  EXPECT_TRUE(contact.contactPoint().isApprox(contact_point_ref));
}


TEST_F(FixedBasePointContactTest, moveAssign) {
  PointContact contact_ref(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  contact_ref.activate();
  const Eigen::Vector3d contact_point_ref = Eigen::Vector3d::Random();
  contact_ref.setContactPoint(contact_point_ref);
  PointContact contact;
  contact = std::move(contact_ref);
  // Check the constructor works well.
  EXPECT_TRUE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), friction_coeff_);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), restitution_coeff_);
  EXPECT_TRUE(contact.contactPoint().isApprox(contact_point_ref));
}


TEST_F(FixedBasePointContactTest, moveConstructor) {
  PointContact contact_ref(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  contact_ref.activate();
  const Eigen::Vector3d contact_point_ref = Eigen::Vector3d::Random();
  contact_ref.setContactPoint(contact_point_ref);
  PointContact contact(std::move(contact_ref));
  // Check the constructor works well.
  EXPECT_TRUE(contact.isActive());
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_DOUBLE_EQ(contact.frictionCoefficient(), friction_coeff_);
  EXPECT_DOUBLE_EQ(contact.restitutionCoefficient(), restitution_coeff_);
  EXPECT_TRUE(contact.contactPoint().isApprox(contact_point_ref));
}


TEST_F(FixedBasePointContactTest, stdvector) {
  std::vector<PointContact> contacts;
  // Check the size at default.
  EXPECT_TRUE(contacts.empty());
  EXPECT_EQ(contacts.size(), 0);
  contacts.push_back(PointContact(model_, contact_frame_id_));
  // Check the size when an element is appended.
  EXPECT_FALSE(contacts.empty());
  EXPECT_EQ(contacts.size(), 1);
  contacts.push_back(PointContact(model_, (int)(contact_frame_id_/2)));
  // Check each contact's parameters.
  EXPECT_EQ(contacts.size(), 2);
  EXPECT_EQ(contacts[0].contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contacts[1].contact_frame_id(), (int)(contact_frame_id_/2));
  int remove_idx = 0;
  std::swap<PointContact>(contacts[remove_idx], contacts.back());
  contacts.pop_back();
  // Check whether swap-and-pop_back works correctly.
  EXPECT_EQ(contacts.size(), 1);
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  contacts.push_back(PointContact(model_, contact_frame_id_));
  EXPECT_EQ(contacts.size(), 2);
  // Check each contact's parameters.
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  EXPECT_EQ(contacts[1].contact_frame_id(), contact_frame_id_);
  remove_idx = 1;
  contacts.erase(contacts.begin()+remove_idx);
  // Check whether erase works correctly.
  EXPECT_EQ(contacts.size(), 1);
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  contacts.push_back(PointContact(model_, contact_frame_id_));
  contacts.push_back(PointContact(model_, contact_frame_id_));
  // Check whether erase works correctly.
  EXPECT_EQ(contacts.size(), 3);
  contacts.erase(contacts.begin()+remove_idx, contacts.end());
  EXPECT_EQ(contacts.size(), 1);
  contacts.pop_back();
  EXPECT_TRUE(contacts.empty());
}


TEST_F(FixedBasePointContactTest, computeJointForceFromContactForce) {
  PointContact contact(model_, contact_frame_id_);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(
            model_.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  Eigen::Vector3d fext = Eigen::Vector3d::Random();
  contact.computeJointForceFromContactForce(fext, fjoint);
  const int parent_joint_id = contact.parent_joint_id();
  fjoint_ref[parent_joint_id] 
      = model_.frames[contact_frame_id_].placement.act(
          pinocchio::Force(fext, Eigen::Vector3d::Zero()));
  for (int i=0; i<fjoint.size(); ++i) {
    EXPECT_TRUE(fjoint[i].isApprox(fjoint_ref[i]));
  }
}


TEST_F(FixedBasePointContactTest, getContactJacobian) {
  PointContact contact(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  pinocchio::forwardKinematics(model_, data_, q_);
  pinocchio::computeJointJacobians(model_, data_, q_);
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd J_trans = Eigen::MatrixXd::Zero(dimv_, 3);
  Eigen::MatrixXd J_ref = Eigen::MatrixXd::Zero(6, dimv_);
  contact.getContactJacobian(model_, data_, J);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL, J_ref);
  EXPECT_TRUE(J.isApprox(J_ref.topRows(3)));
  const bool transpose = true;
  contact.getContactJacobian(model_, data_, J_trans, transpose);
  EXPECT_TRUE(J_trans.isApprox(J_ref.topRows(3).transpose()));
  const int block_rows_begin = rand() % 5;
  const int block_cols_begin = rand() % 5;
  Eigen::MatrixXd J_block 
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  contact.getContactJacobian(model_, data_, J_block.block(block_rows_begin, 
                                                          block_cols_begin,
                                                          3, dimv_));
  EXPECT_TRUE(J_block.block(block_rows_begin, block_cols_begin, 3, dimv_).
              isApprox(J));
  Eigen::MatrixXd J_trans_block 
      = Eigen::MatrixXd::Zero(dimv_+2*block_cols_begin, 3+2*block_rows_begin);
  contact.getContactJacobian(model_, data_, J_trans_block.block(block_cols_begin, 
                                                                block_rows_begin,
                                                                dimv_, 3), true);
  EXPECT_TRUE(J_block.isApprox(J_trans_block.transpose()));
  J_block.block(block_rows_begin, block_cols_begin, 3, dimv_) -= J;
  EXPECT_TRUE(J_block.isZero());
}


TEST_F(FixedBasePointContactTest, baumgarteResidual) {
  PointContact contact(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::Vector3d residual, residual_ref;
  residual.setZero();
  residual_ref.setZero();
  contact.setContactPointByCurrentKinematics(data_);
  const double time_step = 0.5;
  contact.computeBaumgarteResidual(model_, data_, time_step, residual);
  const double baumgarte_weight_on_velocity = (2-restitution_coeff_) / time_step;
  const double baumgarte_weight_on_position = 1 / (time_step*time_step);
  residual_ref 
      = pinocchio::getFrameClassicalAcceleration(model_, data_, contact_frame_id_, 
                                                 pinocchio::LOCAL).linear()
          + baumgarte_weight_on_velocity
              * pinocchio::getFrameVelocity(model_, data_, contact_frame_id_, 
                                            pinocchio::LOCAL).linear()
          + baumgarte_weight_on_position 
              * (data_.oMf[contact_frame_id_].translation()
                 -contact.contactPoint());
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::VectorXd residuals = Eigen::VectorXd::Zero(10);
  contact.computeBaumgarteResidual(model_, data_, time_step, residuals.segment<3>(5));
  EXPECT_TRUE(residuals.head(5).isZero());
  EXPECT_TRUE(residuals.segment<3>(5).isApprox(residual_ref));
  EXPECT_TRUE(residuals.tail(2).isZero());
  const double coeff = Eigen::Vector2d::Random()[0];
  contact.computeBaumgarteResidual(model_, data_, coeff, time_step, residuals.segment<3>(5));
  EXPECT_TRUE(residuals.head(5).isZero());
  EXPECT_TRUE(residuals.segment<3>(5).isApprox(coeff*residual_ref));
  EXPECT_TRUE(residuals.tail(2).isZero());
}


TEST_F(FixedBasePointContactTest, baumgarteDerivatives) {
  PointContact contact(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(3, dimv_);
  const double time_step = 0.5;
  contact.computeBaumgarteDerivatives(model_, data_, time_step, baum_partial_dq, 
                                      baum_partial_dv, baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd frame_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd frame_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd frame_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd frame_a_partial_da = Eigen::MatrixXd::Zero(6, dimv_);
  pinocchio::getFrameAccelerationDerivatives(model_, data_, contact_frame_id_, 
                                             pinocchio::LOCAL,
                                             frame_v_partial_dq, 
                                             frame_a_partial_dq, 
                                             frame_a_partial_dv, 
                                             frame_a_partial_da);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv_);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL, J_frame);
  pinocchio::Motion v_frame = pinocchio::getFrameVelocity(model_, data_, 
                                                          contact_frame_id_, 
                                                          pinocchio::LOCAL);
  Eigen::Matrix3d v_linear_skew, v_angular_skew;
  v_linear_skew.setZero();
  v_angular_skew.setZero();
  pinocchio::skew(v_frame.linear(), v_linear_skew);
  pinocchio::skew(v_frame.angular(), v_angular_skew);
  const double baumgarte_weight_on_velocity = (2-restitution_coeff_) / time_step;
  const double baumgarte_weight_on_position = 1 / (time_step*time_step);
  baum_partial_dq_ref 
      = frame_a_partial_dq.template topRows<3>()
          + v_angular_skew * frame_v_partial_dq.template topRows<3>()
          + v_linear_skew * frame_v_partial_dq.template bottomRows<3>()
          + baumgarte_weight_on_velocity 
              * frame_v_partial_dq.template topRows<3>()
          + baumgarte_weight_on_position 
              * data_.oMf[contact_frame_id_].rotation()
              * J_frame.template topRows<3>();
  baum_partial_dv_ref 
      = frame_a_partial_dv.template topRows<3>()
          + v_angular_skew * J_frame.template topRows<3>()
          + v_linear_skew * J_frame.template bottomRows<3>()
          + baumgarte_weight_on_velocity
              * frame_a_partial_da.template topRows<3>();
  baum_partial_da_ref 
      = frame_a_partial_da.template topRows<3>();
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
  std::cout << "baum_partial_dq:" << std::endl;
  std::cout << baum_partial_dq << std::endl;
  std::cout << std::endl;
  std::cout << "baum_partial_dv:" << std::endl;
  std::cout << baum_partial_dv << std::endl;
  std::cout << std::endl;
  std::cout << "baum_partial_da:" << std::endl;
  std::cout << baum_partial_da << std::endl;
  std::cout << std::endl;
}


TEST_F(FixedBasePointContactTest, baumgarteDerivativesBlock) {
  PointContact contact(model_, contact_frame_id_, friction_coeff_, restitution_coeff_);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(3, dimv_);
  const double time_step = 0.5;
  contact.computeBaumgarteDerivatives(model_, data_, time_step, baum_partial_dq_ref, 
                                      baum_partial_dv_ref, baum_partial_da_ref);
  const int block_rows_begin = rand() % 5;
  const int block_cols_begin = rand() % 5;
  Eigen::MatrixXd baum_partial_dq 
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  Eigen::MatrixXd baum_partial_dv 
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  Eigen::MatrixXd baum_partial_da 
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  contact.computeBaumgarteDerivatives(
      model_, data_, time_step, 
      baum_partial_dq.block(block_rows_begin, block_cols_begin, 3, dimv_),
      baum_partial_dv.block(block_rows_begin, block_cols_begin, 3, dimv_),
      baum_partial_da.block(block_rows_begin, block_cols_begin, 3, dimv_));
  EXPECT_TRUE(
        baum_partial_dq.block(block_rows_begin, block_cols_begin, 3, dimv_)
        .isApprox(baum_partial_dq_ref));
  EXPECT_TRUE(
        baum_partial_dv.block(block_rows_begin, block_cols_begin, 3, dimv_)
        .isApprox(baum_partial_dv_ref));
  EXPECT_TRUE(
        baum_partial_da.block(block_rows_begin, block_cols_begin, 3, dimv_)
        .isApprox(baum_partial_da_ref));
  Eigen::MatrixXd baum_partial_dq_coeff
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  Eigen::MatrixXd baum_partial_dv_coeff
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  Eigen::MatrixXd baum_partial_da_coeff
      = Eigen::MatrixXd::Zero(3+2*block_rows_begin, dimv_+2*block_cols_begin);
  const double coeff = Eigen::VectorXd::Random(1)[0];
  contact.computeBaumgarteDerivatives(
      model_, data_, coeff, time_step,
      baum_partial_dq.block(block_rows_begin, block_cols_begin, 3, dimv_),
      baum_partial_dv.block(block_rows_begin, block_cols_begin, 3, dimv_),
      baum_partial_da.block(block_rows_begin, block_cols_begin, 3, dimv_));
  EXPECT_TRUE(baum_partial_dq_coeff.isApprox(coeff*baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_coeff.isApprox(coeff*baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_coeff.isApprox(coeff*baum_partial_da));
  baum_partial_dq.block(block_rows_begin, block_cols_begin, 3, dimv_) 
      -= baum_partial_dq_ref;
  baum_partial_dv.block(block_rows_begin, block_cols_begin, 3, dimv_) 
      -= baum_partial_dv_ref;
  baum_partial_da.block(block_rows_begin, block_cols_begin, 3, dimv_) 
      -= baum_partial_da_ref;
  EXPECT_TRUE(baum_partial_dq.isZero());
  EXPECT_TRUE(baum_partial_dv.isZero());
  EXPECT_TRUE(baum_partial_da.isZero());
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}