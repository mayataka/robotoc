#include <vector>
#include <string>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "robot/point_contact.hpp"
#include "robot/robot.hpp"


namespace idocp {

class FloatingBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf_ = "../../../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(urdf_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    contact_frames_.clear();
    contact_frames_ = {14, 24, 34, 44};
    q_ = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                        Eigen::VectorXd::Ones(dimq_));
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    baumgarte_weight_on_velocity_ = std::abs(Eigen::VectorXd::Random(2)[0]);
    baumgarte_weight_on_position_ = std::abs(Eigen::VectorXd::Random(2)[0]);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, max_dimf_;
  std::vector<int> contact_frames_;
  Eigen::VectorXd q_, v_, a_;
  double baumgarte_weight_on_velocity_, baumgarte_weight_on_position_;
};


TEST_F(FloatingBaseRobotTest, constructor) {
  // Default constructor
  Robot robot_empty;
  EXPECT_EQ(robot_empty.dimq(), 0);
  EXPECT_EQ(robot_empty.dimv(), 0);
  EXPECT_EQ(robot_empty.dimf(), 0);
  EXPECT_EQ(robot_empty.max_dimf(), 0);
  EXPECT_EQ(robot_empty.dim_passive(), 0);
  EXPECT_EQ(robot_empty.max_point_contacts(), 0);
  EXPECT_FALSE(robot_empty.has_floating_base());
  EXPECT_TRUE(robot_empty.passive_joint_indices().empty());
  Robot robot(urdf_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimv_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.max_dimf(), 0);
  EXPECT_EQ(robot.dim_passive(), 6);
  EXPECT_EQ(robot.max_point_contacts(), 0);
  EXPECT_TRUE(robot.has_floating_base());
  EXPECT_FALSE(robot.passive_joint_indices().empty());
  robot.printRobotModel();
  Robot robot_contact(urdf_, contact_frames_, 
                      baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  EXPECT_EQ(robot_contact.dimq(), dimq_);
  EXPECT_EQ(robot_contact.dimv(), dimv_);
  EXPECT_EQ(robot_contact.dimf(), 0);
  EXPECT_EQ(robot_contact.max_dimf(), 3*contact_frames_.size());
  EXPECT_EQ(robot_contact.dim_passive(), 6);
  EXPECT_EQ(robot_contact.max_point_contacts(), contact_frames_.size());
  EXPECT_TRUE(robot_contact.has_floating_base());
  EXPECT_FALSE(robot_contact.passive_joint_indices().empty());
  robot_contact.printRobotModel();
  Eigen::VectorXd effort_limit, velocity_limit, lower_position_limit, 
                  upper_position_limit;
  effort_limit = Eigen::VectorXd::Constant(dimv_-6, 80);
  velocity_limit = Eigen::VectorXd::Constant(dimv_-6, 15);
  lower_position_limit = Eigen::VectorXd::Constant(dimv_-6, -9.42);
  upper_position_limit = Eigen::VectorXd::Constant(dimv_-6, 9.42);
  EXPECT_TRUE(
      robot.jointEffortLimit().tail(dimv_-6).isApprox(effort_limit));
  EXPECT_TRUE(
      robot.jointVelocityLimit().tail(dimv_-6).isApprox(velocity_limit));
  EXPECT_TRUE(
      robot.lowerJointPositionLimit().tail(dimv_-6)
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot.upperJointPositionLimit().tail(dimv_-6)
      .isApprox(upper_position_limit));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().tail(dimv_-6).isApprox(effort_limit));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit().tail(dimv_-6)
      .isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit().tail(dimv_-6)
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit().tail(dimv_-6)
      .isApprox(upper_position_limit));
}


TEST_F(FloatingBaseRobotTest, integrateConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = q_;
  Eigen::VectorXd q_ref = q_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q);
  q_ref = pinocchio::integrate(model_, q_, integration_length*v_);
  EXPECT_TRUE(q.isApprox(q_ref));
}


TEST_F(FloatingBaseRobotTest, differenceConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = q_;
  Eigen::VectorXd v_ref = v_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q);
  robot.differenceConfiguration(q, q_, v_ref);
  v_ref = v_ref / integration_length;
  EXPECT_TRUE(v_.isApprox(v_ref));
}


TEST_F(FloatingBaseRobotTest, dIntegrateConfiguration) {
  Robot robot(urdf_);
  Eigen::MatrixXd dintegrate_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dintegrate_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dintegrate_dq_ref = dintegrate_dq;
  Eigen::MatrixXd dintegrate_dv_ref = dintegrate_dv;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.dIntegrateConfiguration(q_, v_, integration_length, dintegrate_dq, 
                                dintegrate_dv);
  pinocchio::dIntegrate(model_, q_, integration_length*v_, dintegrate_dq_ref, 
                        pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q_, integration_length*v_, dintegrate_dv_ref, 
                        pinocchio::ARG1);
  EXPECT_TRUE(dintegrate_dq.isApprox(dintegrate_dq_ref));
  EXPECT_TRUE(dintegrate_dv.isApprox(dintegrate_dv_ref));
  EXPECT_TRUE(
      dintegrate_dq.block(6, 6, dimv_-6, dimv_-6)
      .isApprox(Eigen::MatrixXd::Identity(dimv_-6, dimv_-6)));
  EXPECT_TRUE(
      dintegrate_dv.block(6, 6, dimv_-6, dimv_-6)
      .isApprox(Eigen::MatrixXd::Identity(dimv_-6, dimv_-6)));
  std::cout << "dintegrate_dq:" << std::endl;
  std::cout << dintegrate_dq << std::endl;
  std::cout << std::endl;
  std::cout << "dintegrate_dv:" << std::endl;
  std::cout << dintegrate_dv << std::endl;
  std::cout << std::endl;
}


TEST_F(FloatingBaseRobotTest, baumgarteResidualAndDerivatives) {
  std::vector<PointContact> contacts_ref; 
  for (const auto& frame : contact_frames_) {
    contacts_ref.push_back(PointContact(model_, frame, 
                                        baumgarte_weight_on_velocity_,
                                        baumgarte_weight_on_position_));
  }
  Robot robot(urdf_, contact_frames_, baumgarte_weight_on_velocity_, 
              baumgarte_weight_on_position_);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  std::vector<bool> is_each_contacts_active(contacts_ref.size(), true);
  robot.setActiveContacts(is_each_contacts_active);
  EXPECT_EQ(robot.dimf(), robot.max_dimf());
  robot.updateKinematics(q_, v_, a_);
  robot.computeBaumgarteResidual(residual);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].resetContactPointByCurrentKinematics(data_);
  }
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeBaumgarteResidual(model_, data_, 3*i, residual_ref);
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));
  std::random_device rnd;
  const int block_rows_begin = rnd() % 10;
  Eigen::MatrixXd baumgarte_partial_q 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimv_);
  Eigen::MatrixXd baumgarte_partial_v 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimv_);
  Eigen::MatrixXd baumgarte_partial_a 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimv_);
  Eigen::MatrixXd baumgarte_partial_q_ref
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimv_);
  Eigen::MatrixXd baumgarte_partial_v_ref 
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimv_);
  Eigen::MatrixXd baumgarte_partial_a_ref 
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimv_);
  robot.computeBaumgarteDerivatives(block_rows_begin, baumgarte_partial_q, 
                                    baumgarte_partial_v, baumgarte_partial_a);
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeBaumgarteDerivatives(model_, data_, 3*i,
                                                baumgarte_partial_q_ref, 
                                                baumgarte_partial_v_ref, 
                                                baumgarte_partial_a_ref);
  }
  EXPECT_TRUE(baumgarte_partial_q.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(baumgarte_partial_v.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(baumgarte_partial_a.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(
      baumgarte_partial_q.bottomRows(robot.max_dimf())
      .isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(
      baumgarte_partial_v.bottomRows(robot.max_dimf())
      .isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(
      baumgarte_partial_a.bottomRows(robot.max_dimf())
      .isApprox(baumgarte_partial_a_ref));
}


TEST_F(FloatingBaseRobotTest, RNEA) {
  // without contact
  Robot robot(urdf_);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(dimv_);
  robot.RNEA(q_, v_, a_, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  Robot robot_contact(urdf_, contact_frames_, baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  tau = Eigen::VectorXd::Zero(dimv_);
  tau_ref = Eigen::VectorXd::Zero(dimv_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot_contact.max_dimf());
  robot_contact.RNEA(q_, v_, a_, tau);
  tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  // with contact
  std::vector<PointContact> contacts_ref; 
  for (const auto& frame : contact_frames_) {
    contacts_ref.push_back(PointContact(model_, frame, 
                                        baumgarte_weight_on_velocity_,
                                        baumgarte_weight_on_position_));
  }
  std::vector<bool> is_each_contacts_active(contacts_ref.size(), true);
  robot_contact.setActiveContacts(is_each_contacts_active);
  robot_contact.setContactForces(fext);
  robot_contact.RNEA(q_, v_, a_, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeJointForceFromContactForce(fext.segment<3>(3*i), 
                                                      fjoint);
  }
  tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));
}


TEST_F(FloatingBaseRobotTest, RNEADerivativesWithoutFext) {
  Robot robot(urdf_);
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = dRNEA_da;
  robot.RNEADerivatives(q_, v_, a_, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::computeRNEADerivatives(model_, data_, q_, v_, a_, dRNEA_dq_ref, 
                                    dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
}


TEST_F(FloatingBaseRobotTest, RNEADerivativesWithContacts) {
  Robot robot(urdf_, contact_frames_, baumgarte_weight_on_velocity_, 
              baumgarte_weight_on_position_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot.max_dimf());
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_dfext = Eigen::MatrixXd::Zero(dimv_, robot.max_dimf());
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = Eigen::MatrixXd::Zero(dimv_, dimv_);
  Eigen::MatrixXd dRNEA_dfext_ref 
      = Eigen::MatrixXd::Zero(dimv_, robot.max_dimf());
  std::vector<PointContact> contacts_ref; 
  for (const auto& frame : contact_frames_) {
    contacts_ref.push_back(PointContact(model_, frame, 
                                        baumgarte_weight_on_velocity_,
                                        baumgarte_weight_on_position_));
  }
  std::vector<bool> is_each_contacts_active(contacts_ref.size(), true);
  robot.setActiveContacts(is_each_contacts_active);
  robot.setContactForces(fext);
  robot.RNEADerivatives(q_, v_, a_, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeJointForceFromContactForce(fext.segment<3>(3*i), 
                                                       fjoint);
  }
  pinocchio::computeRNEADerivatives(model_, data_, q_, v_, a_, fjoint, 
                                    dRNEA_dq_ref, dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
  robot.dRNEAPartialdFext(dRNEA_dfext);
  const bool transpose_jacobian = true;
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].getContactJacobian(model_, data_, 3*i, dRNEA_dfext_ref, 
                                       transpose_jacobian);
  }
  EXPECT_TRUE(dRNEA_dfext.isApprox(dRNEA_dfext_ref));
}


TEST_F(FloatingBaseRobotTest, floating_base) {
  Robot robot(urdf_);
  Eigen::VectorXd tau = Eigen::VectorXd::Ones(robot.dimv());
  robot.setPassiveTorques(tau);
  EXPECT_TRUE(tau.head(6).isApprox(Eigen::VectorXd::Zero(6)));
  EXPECT_TRUE(tau.tail(dimv_-6).isApprox(Eigen::VectorXd::Ones(dimv_-6)));
}

} // namespace idocp 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}