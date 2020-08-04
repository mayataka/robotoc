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

class FixedBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf_ = "../../urdf/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    contact_frame_id_ = 18;
    q_ = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                        Eigen::VectorXd::Ones(dimq_));
    v_ = Eigen::VectorXd::Random(dimq_);
    a_ = Eigen::VectorXd::Random(dimq_);
    baumgarte_weight_on_velocity_ = std::abs(Eigen::VectorXd::Random(2)[0]);
    baumgarte_weight_on_position_ = std::abs(Eigen::VectorXd::Random(2)[0]);
  }

  virtual void TearDown() {
  }

  std::string urdf_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, contact_frame_id_;
  Eigen::VectorXd q_, v_, a_;
  double baumgarte_weight_on_velocity_, baumgarte_weight_on_position_;
};


TEST_F(FixedBaseRobotTest, constructor) {
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
  EXPECT_EQ(robot.dimv(), dimq_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.max_dimf(), 0);
  EXPECT_EQ(robot.dim_passive(), 0);
  EXPECT_EQ(robot.max_point_contacts(), 0);
  EXPECT_FALSE(robot.has_floating_base());
  robot.printRobotModel();
  std::vector<int> contacts = {contact_frame_id_};
  Robot robot_contact(urdf_, contacts, baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  EXPECT_EQ(robot_contact.dimq(), dimq_);
  EXPECT_EQ(robot_contact.dimv(), dimq_);
  EXPECT_EQ(robot_contact.dimf(), 0);
  EXPECT_EQ(robot_contact.max_dimf(), 3);
  EXPECT_EQ(robot_contact.dim_passive(), 0);
  EXPECT_EQ(robot_contact.max_point_contacts(), 1);
  EXPECT_EQ(robot_contact.is_contact_active(0), false);
  EXPECT_FALSE(robot_contact.has_floating_base());
  robot_contact.printRobotModel();
  Eigen::VectorXd effort_limit(dimq_), velocity_limit(dimq_), 
                  lower_position_limit(dimq_), upper_position_limit(dimq_);
  effort_limit << 300, 300, 300, 300, 300, 300, 300;
  velocity_limit << 10, 10, 10, 10, 10, 10, 10;
  lower_position_limit << -2.96705972839, -2.09439510239, -2.96705972839, 
                          -2.09439510239, -2.96705972839, -2.09439510239, 
                          -3.05432619099;
  upper_position_limit << 2.96705972839, 2.09439510239, 2.96705972839, 
                          2.09439510239, 2.96705972839, 2.09439510239, 
                          3.05432619099;
  EXPECT_TRUE(robot.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(robot.lowerJointPositionLimit().isApprox(lower_position_limit));
  EXPECT_TRUE(robot.upperJointPositionLimit().isApprox(upper_position_limit));
  EXPECT_TRUE(robot_contact.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot_contact.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(
      robot_contact.lowerJointPositionLimit().isApprox(lower_position_limit));
  EXPECT_TRUE(
      robot_contact.upperJointPositionLimit().isApprox(upper_position_limit));
  effort_limit.setZero();
  velocity_limit.setZero();
  lower_position_limit.setZero();
  upper_position_limit.setZero();
  robot.setJointEffortLimit(effort_limit);
  robot.setJointVelocityLimit(velocity_limit);
  robot.setLowerJointPositionLimit(lower_position_limit);
  robot.setUpperJointPositionLimit(upper_position_limit);
  EXPECT_TRUE(robot.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(robot.lowerJointPositionLimit().isApprox(lower_position_limit));
  EXPECT_TRUE(robot.upperJointPositionLimit().isApprox(upper_position_limit));
}


TEST_F(FixedBaseRobotTest, integrateConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = q_;
  Eigen::VectorXd q_ref = q_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q);
  q_ref = pinocchio::integrate(model_, q_, integration_length*v_);
  EXPECT_TRUE(q.isApprox(q_ref));
}


TEST_F(FixedBaseRobotTest, differenceConfiguration) {
  Robot robot(urdf_);
  Eigen::VectorXd q = q_;
  Eigen::VectorXd v_ref = v_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q);
  robot.differenceConfiguration(q, q_, v_ref);
  v_ref = v_ref / integration_length;
  EXPECT_TRUE(v_.isApprox(v_ref));
}


TEST_F(FixedBaseRobotTest, dIntegrateConfiguration) {
  Robot robot(urdf_);
  Eigen::MatrixXd dintegrate_dq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dintegrate_dv = Eigen::MatrixXd::Zero(dimq_, dimq_);
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.dIntegrateConfiguration(q_, v_, integration_length, dintegrate_dq, 
                                dintegrate_dv);
  EXPECT_TRUE(dintegrate_dq.isApprox(Eigen::MatrixXd::Identity(dimq_, dimq_)));
  EXPECT_TRUE(dintegrate_dv.isApprox(Eigen::MatrixXd::Identity(dimq_, dimq_)));
}


TEST_F(FixedBaseRobotTest, configurationJacobian) {
  Robot robot(urdf_);
  Eigen::MatrixXd Jacobian = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd Jacobian_ref = Eigen::MatrixXd::Zero(dimq_, dimq_);
  robot.configurationJacobian(q_, Jacobian);
  pinocchio::integrateCoeffWiseJacobian(model_, q_, Jacobian_ref);
  EXPECT_TRUE(Jacobian.isApprox(Jacobian_ref));
  EXPECT_TRUE(Jacobian.isApprox(Eigen::MatrixXd::Identity(dimq_, dimq_)));
  std::cout << "configuration Jacobian:" << std::endl;
  std::cout << Jacobian << std::endl;
  std::cout << std::endl;
}


TEST_F(FixedBaseRobotTest, baumgarteResidualAndDerivatives) {
  std::vector<int> contact_frames = {contact_frame_id_};
  Robot robot(urdf_, contact_frames, baumgarte_weight_on_velocity_, 
              baumgarte_weight_on_position_);
  std::random_device rnd;
  const int block_begin = rnd() % 10;
  Eigen::VectorXd residual 
      = Eigen::VectorXd::Zero(block_begin+robot.max_dimf());
  Eigen::VectorXd residual_ref 
      = Eigen::VectorXd::Zero(block_begin+robot.max_dimf());
  std::vector<bool> is_each_contacts_active = {true};
  robot.setActiveContacts(is_each_contacts_active);
  EXPECT_EQ(robot.dimf(), robot.max_dimf());
  EXPECT_EQ(robot.is_contact_active(0), true);
  robot.updateKinematics(q_, v_, a_);
  robot.computeBaumgarteResidual(block_begin, residual);
  PointContact contact_ref(model_, contact_frame_id_, 
                           baumgarte_weight_on_velocity_, 
                           baumgarte_weight_on_position_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  contact_ref.resetContactPointByCurrentKinematics(data_);
  contact_ref.computeBaumgarteResidual(model_, data_, block_begin, residual_ref);
  EXPECT_TRUE(residual.isApprox(residual_ref));
  const double coeff = Eigen::VectorXd::Random(1)[0];
  robot.computeBaumgarteResidual(block_begin, coeff, residual);
  EXPECT_TRUE(residual.isApprox(coeff*residual_ref));
  const int block_rows_begin = rnd() % 10;
  Eigen::MatrixXd baumgarte_partial_q 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_v 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_a 
      = Eigen::MatrixXd::Zero(block_rows_begin+robot.max_dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_q_ref
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_v_ref 
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_a_ref 
      = Eigen::MatrixXd::Zero(robot.max_dimf(), dimq_);
  robot.computeBaumgarteDerivatives(block_rows_begin, baumgarte_partial_q, 
                                    baumgarte_partial_v, baumgarte_partial_a);
  contact_ref.computeBaumgarteDerivatives(model_, data_, 
                                          baumgarte_partial_q_ref, 
                                          baumgarte_partial_v_ref, 
                                          baumgarte_partial_a_ref);
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
  baumgarte_partial_q.setZero();
  baumgarte_partial_v.setZero();
  baumgarte_partial_a.setZero();
  robot.computeBaumgarteDerivatives(block_rows_begin, coeff, baumgarte_partial_q, 
                                    baumgarte_partial_v, baumgarte_partial_a);
  EXPECT_TRUE(baumgarte_partial_q.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(baumgarte_partial_v.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(baumgarte_partial_a.topRows(block_rows_begin).isZero());
  EXPECT_TRUE(
      baumgarte_partial_q.bottomRows(robot.max_dimf())
      .isApprox(coeff*baumgarte_partial_q_ref));
  EXPECT_TRUE(
      baumgarte_partial_v.bottomRows(robot.max_dimf())
      .isApprox(coeff*baumgarte_partial_v_ref));
  EXPECT_TRUE(
      baumgarte_partial_a.bottomRows(robot.max_dimf())
      .isApprox(coeff*baumgarte_partial_a_ref));
}


TEST_F(FixedBaseRobotTest, RNEA) {
  // without contact
  Robot robot(urdf_);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(dimq_);
  robot.RNEA(q_, v_, a_, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  std::vector<int> contact_frames = {contact_frame_id_};
  Robot robot_contact(urdf_, contact_frames, baumgarte_weight_on_velocity_, 
                      baumgarte_weight_on_position_);
  tau = Eigen::VectorXd::Zero(dimq_);
  tau_ref = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot_contact.max_dimf());
  robot_contact.RNEA(q_, v_, a_, tau);
  tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  // with contact
  std::vector<bool> is_each_contacts_active = {true};
  robot_contact.setActiveContacts(is_each_contacts_active);
  robot_contact.setContactForces(fext);
  robot_contact.RNEA(q_, v_, a_, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, 
                           baumgarte_weight_on_velocity_, 
                           baumgarte_weight_on_position_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));
}


TEST_F(FixedBaseRobotTest, RNEADerivativesWithoutFext) {
  Robot robot(urdf_);
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(dimq_, dimq_);
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


TEST_F(FixedBaseRobotTest, RNEADerivativesWithContacts) {
  std::vector<int> contact_frames = {contact_frame_id_};
  Robot robot(urdf_, contact_frames, baumgarte_weight_on_velocity_, 
              baumgarte_weight_on_position_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot.max_dimf());
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dfext = Eigen::MatrixXd::Zero(dimq_, robot.max_dimf());
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dfext_ref 
      = Eigen::MatrixXd::Zero(dimq_, robot.max_dimf());
  std::vector<bool> is_each_contacts_active = {true};
  robot.setActiveContacts(is_each_contacts_active);
  robot.setContactForces(fext);
  robot.RNEADerivatives(q_, v_, a_, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, 
                           baumgarte_weight_on_velocity_, 
                           baumgarte_weight_on_position_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  pinocchio::computeRNEADerivatives(model_, data_, q_, v_, a_, fjoint, 
                                    dRNEA_dq_ref, dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
  const bool transpose_jacobian = true;
  robot.dRNEAPartialdFext(dRNEA_dfext);
  contact_ref.getContactJacobian(model_, data_, dRNEA_dfext_ref, 
                                 transpose_jacobian);
  EXPECT_TRUE(dRNEA_dfext.isApprox(dRNEA_dfext_ref));
}


TEST_F(FixedBaseRobotTest, floating_base) {
  Robot robot(urdf_);
  Eigen::VectorXd tau = Eigen::VectorXd::Ones(robot.dimv());
  robot.setPassiveTorques(tau);
  EXPECT_TRUE(tau.isApprox(Eigen::VectorXd::Ones(robot.dimv())));
}

} // namespace idocp 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}