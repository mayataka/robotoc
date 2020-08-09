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

#include "idocp/robot/point_contact.hpp"
#include "idocp/robot/robot.hpp"


namespace idocp {

class FloatingBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf_ = "../../urdf/anymal/anymal.urdf";
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
  EXPECT_EQ(robot_empty.num_active_point_contacts(), 0);
  EXPECT_FALSE(robot_empty.has_floating_base());
  Robot robot(urdf_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimv_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.max_dimf(), 0);
  EXPECT_EQ(robot.dim_passive(), 6);
  EXPECT_EQ(robot.max_point_contacts(), 0);
  EXPECT_EQ(robot.num_active_point_contacts(), 0);
  EXPECT_TRUE(robot.has_floating_base());
  robot.printRobotModel();
  Robot robot_contact(urdf_, contact_frames_, 
                      baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  EXPECT_EQ(robot_contact.dimq(), dimq_);
  EXPECT_EQ(robot_contact.dimv(), dimv_);
  EXPECT_EQ(robot_contact.dimf(), 0);
  EXPECT_EQ(robot_contact.max_dimf(), 3*contact_frames_.size());
  EXPECT_EQ(robot_contact.dim_passive(), 6);
  for (int i=0; i<robot_contact.max_point_contacts(); ++i) {
    EXPECT_FALSE(robot_contact.is_contact_active(i));
  }
  EXPECT_EQ(robot_contact.max_point_contacts(), contact_frames_.size());
  EXPECT_EQ(robot_contact.num_active_point_contacts(), 0);
  EXPECT_TRUE(robot_contact.has_floating_base());
  robot_contact.printRobotModel();
  Eigen::VectorXd effort_limit, velocity_limit, lower_position_limit, 
                  upper_position_limit;
  effort_limit = Eigen::VectorXd::Constant(dimv_-6, 80);
  velocity_limit = Eigen::VectorXd::Constant(dimv_-6, 15);
  lower_position_limit = Eigen::VectorXd::Constant(dimv_-6, -9.42);
  upper_position_limit = Eigen::VectorXd::Constant(dimv_-6, 9.42);
  EXPECT_TRUE(robot.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(
      robot.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit()
      .isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  robot.setJointEffortLimit(effort_limit);
  robot.setJointVelocityLimit(velocity_limit);
  robot.setLowerJointPositionLimit(lower_position_limit);
  robot.setUpperJointPositionLimit(upper_position_limit);
  EXPECT_TRUE(robot.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(robot.lowerJointPositionLimit().isApprox(lower_position_limit));
  EXPECT_TRUE(robot.upperJointPositionLimit().isApprox(upper_position_limit));
}


TEST_F(FloatingBaseRobotTest, moveAssign) {
  // Default constructor
  Robot robot_empty;
  EXPECT_EQ(robot_empty.dimq(), 0);
  EXPECT_EQ(robot_empty.dimv(), 0);
  EXPECT_EQ(robot_empty.dimf(), 0);
  EXPECT_EQ(robot_empty.max_dimf(), 0);
  EXPECT_EQ(robot_empty.dim_passive(), 0);
  EXPECT_EQ(robot_empty.max_point_contacts(), 0);
  EXPECT_FALSE(robot_empty.has_floating_base());
  Robot robot(urdf_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimv_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.max_dimf(), 0);
  EXPECT_EQ(robot.dim_passive(), 6);
  EXPECT_EQ(robot.max_point_contacts(), 0);
  EXPECT_TRUE(robot.has_floating_base());
  robot.printRobotModel();
  Robot robot_contact(urdf_, contact_frames_, 
                      baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  EXPECT_EQ(robot_contact.dimq(), dimq_);
  EXPECT_EQ(robot_contact.dimv(), dimv_);
  EXPECT_EQ(robot_contact.dimf(), 0);
  EXPECT_EQ(robot_contact.max_dimf(), 3*contact_frames_.size());
  EXPECT_EQ(robot_contact.dim_passive(), 6);
  for (int i=0; i<robot_contact.max_point_contacts(); ++i) {
    EXPECT_EQ(robot_contact.is_contact_active(i), false);
  }
  EXPECT_EQ(robot_contact.max_point_contacts(), contact_frames_.size());
  EXPECT_TRUE(robot_contact.has_floating_base());
  robot_contact.printRobotModel();
  Eigen::VectorXd effort_limit, velocity_limit, lower_position_limit, 
                  upper_position_limit;
  effort_limit = Eigen::VectorXd::Constant(dimv_-6, 80);
  velocity_limit = Eigen::VectorXd::Constant(dimv_-6, 15);
  lower_position_limit = Eigen::VectorXd::Constant(dimv_-6, -9.42);
  upper_position_limit = Eigen::VectorXd::Constant(dimv_-6, 9.42);
  EXPECT_TRUE(robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot_contact.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit()
      .isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  Robot robot_ref = robot_contact; 
  robot_empty = std::move(robot_ref);
  EXPECT_EQ(robot_contact.dimq(), robot_empty.dimq());
  EXPECT_EQ(robot_contact.dimv(), robot_empty.dimv());
  EXPECT_EQ(robot_contact.dimf(), robot_empty.dimf());
  EXPECT_EQ(robot_contact.max_dimf(), robot_empty.max_dimf());
  EXPECT_EQ(robot_contact.dim_passive(), robot_empty.dim_passive());
  EXPECT_EQ(robot_contact.max_point_contacts(), robot_empty.max_point_contacts());
  EXPECT_TRUE(robot_empty.has_floating_base());
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(robot_empty.jointEffortLimit()));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(robot_empty.jointEffortLimit()));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit().isApprox(robot_empty.jointVelocityLimit()));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit().isApprox(robot_empty.lowerJointPositionLimit()));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit().isApprox(robot_empty.upperJointPositionLimit()));
}


TEST_F(FloatingBaseRobotTest, moveConstructor) {
  // Default constructor
  Robot robot_contact(urdf_, contact_frames_, 
                      baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  EXPECT_EQ(robot_contact.dimq(), dimq_);
  EXPECT_EQ(robot_contact.dimv(), dimv_);
  EXPECT_EQ(robot_contact.dimf(), 0);
  EXPECT_EQ(robot_contact.max_dimf(), 3*contact_frames_.size());
  EXPECT_EQ(robot_contact.dim_passive(), 6);
  for (int i=0; i<robot_contact.max_point_contacts(); ++i) {
    EXPECT_EQ(robot_contact.is_contact_active(i), false);
  }
  EXPECT_EQ(robot_contact.max_point_contacts(), contact_frames_.size());
  EXPECT_TRUE(robot_contact.has_floating_base());
  robot_contact.printRobotModel();
  Eigen::VectorXd effort_limit, velocity_limit, lower_position_limit, 
                  upper_position_limit;
  effort_limit = Eigen::VectorXd::Constant(dimv_-6, 80);
  velocity_limit = Eigen::VectorXd::Constant(dimv_-6, 15);
  lower_position_limit = Eigen::VectorXd::Constant(dimv_-6, -9.42);
  upper_position_limit = Eigen::VectorXd::Constant(dimv_-6, 9.42);
  EXPECT_TRUE(robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot_contact.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit()
      .isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit()
      .isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit()
      .isApprox(upper_position_limit));
  Robot robot_ref = robot_contact; 
  Robot robot_empty(std::move(robot_ref));
  EXPECT_EQ(robot_contact.dimq(), robot_empty.dimq());
  EXPECT_EQ(robot_contact.dimv(), robot_empty.dimv());
  EXPECT_EQ(robot_contact.dimf(), robot_empty.dimf());
  EXPECT_EQ(robot_contact.max_dimf(), robot_empty.max_dimf());
  EXPECT_EQ(robot_contact.dim_passive(), robot_empty.dim_passive());
  EXPECT_EQ(robot_contact.max_point_contacts(), robot_empty.max_point_contacts());
  EXPECT_TRUE(robot_empty.has_floating_base());
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(robot_empty.jointEffortLimit()));
  EXPECT_TRUE(
      robot_contact.jointEffortLimit().isApprox(robot_empty.jointEffortLimit()));
  EXPECT_TRUE(
      robot_contact.jointVelocityLimit().isApprox(robot_empty.jointVelocityLimit()));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit().isApprox(robot_empty.lowerJointPositionLimit()));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit().isApprox(robot_empty.upperJointPositionLimit()));
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


TEST_F(FloatingBaseRobotTest, subtractConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = q_;
  Eigen::VectorXd v_ref = v_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q);
  robot.subtractConfiguration(q, q_, v_ref);
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
  pinocchio::dIntegrate(model_, q_, v_, dintegrate_dq_ref, pinocchio::ARG0);
  pinocchio::dIntegrate(model_, q_, v_, dintegrate_dv_ref, pinocchio::ARG1);
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


TEST_F(FloatingBaseRobotTest, computeConfigurationJacobian) {
  Robot robot(urdf_);
  Eigen::MatrixXd Jacobian_ref = Eigen::MatrixXd::Zero(dimq_, dimv_);
  pinocchio::integrateCoeffWiseJacobian(model_, q_, Jacobian_ref);
  Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(dimq_, dimv_);
  robot.computeConfigurationJacobian(q_, Jacobian);
  EXPECT_TRUE(Jacobian.isApprox(Jacobian_ref));
  std::cout << "configuration Jacobian:" << std::endl;
  std::cout << Jacobian << std::endl;
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
  std::random_device rnd;
  const int segment_begin = rnd() % 5;
  Eigen::VectorXd residual 
      = Eigen::VectorXd::Zero(segment_begin+robot.max_dimf());
  Eigen::VectorXd residual_ref 
      = Eigen::VectorXd::Zero(segment_begin+robot.max_dimf());
  std::vector<bool> is_each_contacts_active(contacts_ref.size(), true);
  robot.setContactStatus(is_each_contacts_active);
  EXPECT_EQ(robot.dimf(), robot.max_dimf());
  for (int i=0; i<robot.max_point_contacts(); ++i) {
    EXPECT_EQ(robot.is_contact_active(i), true);
  }
  robot.updateKinematics(q_, v_, a_);
  robot.setContactPointsByCurrentKinematics();
  robot.computeBaumgarteResidual(segment_begin, residual);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].resetContactPointByCurrentKinematics(data_);
  }
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeBaumgarteResidual(model_, data_, segment_begin+3*i, 
                                             residual_ref);
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));
  const double coeff = Eigen::VectorXd::Random(1)[0];
  robot.computeBaumgarteResidual(segment_begin, coeff, residual);
  EXPECT_TRUE(residual.isApprox(coeff*residual_ref));

  std::cout << "aaa" << std::endl;
  const int block_rows_begin = rnd() % 5;
  const int block_cols_begin = rnd() % 5;
  Eigen::MatrixXd baumgarte_partial_q_ref
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_v_ref 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_a_ref 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  for (int i=0; i<contacts_ref.size(); ++i) {
    contacts_ref[i].computeBaumgarteDerivatives(model_, data_, 
                                                block_rows_begin+3*i, 
                                                block_cols_begin, 
                                                baumgarte_partial_q_ref, 
                                                baumgarte_partial_v_ref, 
                                                baumgarte_partial_a_ref);
  }
  Eigen::MatrixXd baumgarte_partial_q 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_v 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_a 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  robot.computeBaumgarteDerivatives(block_rows_begin, block_cols_begin, baumgarte_partial_q, 
                                    baumgarte_partial_v, baumgarte_partial_a);
  EXPECT_TRUE(baumgarte_partial_q.isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(baumgarte_partial_v.isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(baumgarte_partial_a.isApprox(baumgarte_partial_a_ref));
  Eigen::MatrixXd baumgarte_partial_q_coeff
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_v_coeff 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  Eigen::MatrixXd baumgarte_partial_a_coeff 
      = Eigen::MatrixXd::Zero(2*block_rows_begin+robot.max_dimf(), 2*block_cols_begin+dimv_);
  robot.computeBaumgarteDerivatives(block_rows_begin, block_cols_begin, coeff,  
                                    baumgarte_partial_q_coeff, 
                                    baumgarte_partial_v_coeff, 
                                    baumgarte_partial_a_coeff);
  EXPECT_TRUE(baumgarte_partial_q_coeff.isApprox(coeff*baumgarte_partial_q));
  EXPECT_TRUE(baumgarte_partial_v_coeff.isApprox(coeff*baumgarte_partial_v));
  EXPECT_TRUE(baumgarte_partial_a_coeff.isApprox(coeff*baumgarte_partial_a));
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
  robot_contact.setContactStatus(is_each_contacts_active);
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
  robot.setContactStatus(is_each_contacts_active);
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
    contacts_ref[i].getContactJacobian(model_, data_, 0, 3*i, -1, dRNEA_dfext_ref,
                                       transpose_jacobian);
  }
  EXPECT_TRUE(dRNEA_dfext.isApprox(dRNEA_dfext_ref));
}


TEST_F(FloatingBaseRobotTest, generateFeasibleConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  for (int i=0; i<robot.dimq()-robot.dim_passive()-1; ++i) {
    EXPECT_TRUE(q(robot.dim_passive()+1+i) >= qmin(i));
    EXPECT_TRUE(q(robot.dim_passive()+1+i) <= qmax(i));
  }
  q = Eigen::VectorXd::Random(robot.dimq());
  Eigen::VectorXd q_ref = q;
  robot.normalizeConfiguration(q);
  pinocchio::normalize(model_, q_ref);
  EXPECT_TRUE(q.isApprox(q_ref));
}

} // namespace idocp 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}