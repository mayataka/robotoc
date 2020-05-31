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


namespace invdynocp {

class FixedBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    urdf_file_name_ = "../../urdf/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    contact_frame_id_ = 18;
    max_point_contacts_ = 1;
    q_ = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                        Eigen::VectorXd::Ones(dimq_));
    v_ = Eigen::VectorXd::Random(dimq_);
    a_ = Eigen::VectorXd::Random(dimq_);
    baumgarte_alpha_ = std::abs(Eigen::VectorXd::Random(2)[0]);
    baumgarte_beta_ = std::abs(Eigen::VectorXd::Random(2)[0]);
  }

  virtual void TearDown() {
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, contact_frame_id_, max_point_contacts_;
  Eigen::VectorXd q_, v_, a_;
  double baumgarte_alpha_, baumgarte_beta_;
};


TEST_F(FixedBaseRobotTest, params) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimq_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.dimfmax(), 3*max_point_contacts_);
  EXPECT_EQ(robot.dim_passive(), 0);
  EXPECT_EQ(robot.max_point_contacts(), max_point_contacts_);
  EXPECT_EQ(robot.urdf_file_name(), urdf_file_name_);
}


TEST_F(FixedBaseRobotTest, integrateConfiguration) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  Eigen::VectorXd q_ref = q_;
  Eigen::VectorXd q_ref2 = q_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q_);
  q_ref = pinocchio::integrate(model_, q_ref, integration_length*v_);
  EXPECT_TRUE(q_.isApprox(q_ref));
  robot.integrateConfiguration(q_ref2, v_, integration_length, q_);
  EXPECT_TRUE(q_.isApprox(q_ref));
}


TEST_F(FixedBaseRobotTest, differenceConfiguration) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  Eigen::VectorXd q1 = q_;
  Eigen::VectorXd v_ref = v_;
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[0]);
  robot.integrateConfiguration(v_, integration_length, q1);
  robot.differenceConfigurations(q1, q_, v_ref);
  v_ref = v_ref / integration_length;
  EXPECT_TRUE(v_.isApprox(v_ref));
}


TEST_F(FixedBaseRobotTest, baumgarteResidualAndDerivatives) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.dimf());
  robot.updateKinematics(q_, v_, a_);
  robot.computeBaumgarteResidual(residual);
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  contact_ref.resetContactPointByCurrentKinematics(data_);
  contact_ref.computeBaumgarteResidual(model_, data_, 0, residual_ref);
  EXPECT_TRUE(residual.isApprox(residual_ref));

  Eigen::MatrixXd baumgarte_partial_q = Eigen::MatrixXd::Zero(robot.dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_v = Eigen::MatrixXd::Zero(robot.dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_a = Eigen::MatrixXd::Zero(robot.dimf(), dimq_);
  Eigen::MatrixXd baumgarte_partial_q_ref = baumgarte_partial_q;
  Eigen::MatrixXd baumgarte_partial_v_ref = baumgarte_partial_v;
  Eigen::MatrixXd baumgarte_partial_a_ref = baumgarte_partial_a;
  robot.computeBaumgarteDerivatives(baumgarte_partial_q, baumgarte_partial_v, 
                                    baumgarte_partial_a);
  contact_ref.computeBaumgarteDerivatives(model_, data_, baumgarte_partial_q_ref, 
                                          baumgarte_partial_v_ref, 
                                          baumgarte_partial_a_ref);
  EXPECT_TRUE(baumgarte_partial_q.isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(baumgarte_partial_v.isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(baumgarte_partial_a.isApprox(baumgarte_partial_a_ref));
}


TEST_F(FixedBaseRobotTest, RNEAWithoutFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(dimq_);
  robot.RNEA(q_, v_, a_, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_);
  EXPECT_TRUE(tau_ref.isApprox(tau));
}


TEST_F(FixedBaseRobotTest, RNEAWithFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd tau_ref = Eigen::VectorXd::Zero(dimq_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot.dimf());
  robot.RNEA(q_, v_, a_, tau_ref);
  robot.RNEA(q_, v_, a_, fext, tau);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  robot.setFext(fext);
  robot.RNEA(q_, v_, a_, fext, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  tau_ref = pinocchio::rnea(model_, data_, q_, v_, a_, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));
}


TEST_F(FixedBaseRobotTest, RNEADerivativesWithoutFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
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


TEST_F(FixedBaseRobotTest, RNEADerivativesWithFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  Eigen::VectorXd fext = Eigen::VectorXd::Random(robot.dimf());
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_da_and_fext = Eigen::MatrixXd::Zero(dimq_+robot.dimf(), dimq_);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::MatrixXd dRNEA_dfext_ref = Eigen::MatrixXd::Zero(robot.dimf(), dimq_);
  robot.setFext(fext);
  robot.RNEADerivatives(q_, v_, a_, fext, dRNEA_dq, dRNEA_dv, dRNEA_da_and_fext);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  pinocchio::computeRNEADerivatives(model_, data_, q_, v_, a_, fjoint, 
                                    dRNEA_dq_ref, dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  contact_ref.getContactJacobian(model_, data_, dRNEA_dfext_ref);
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da_and_fext.topRows(dimq_).isApprox(dRNEA_da_ref));
  EXPECT_TRUE(dRNEA_da_and_fext.block(dimq_, 0, 3, dimq_).isApprox(dRNEA_dfext_ref));
}


TEST_F(FixedBaseRobotTest, passiveJoints) {
  Robot robot(urdf_file_name_, 0);
  Eigen::VectorXd tau = Eigen::VectorXd::Ones(robot.dimv());
  robot.setPassiveTorques(tau);
  EXPECT_TRUE(tau.isApprox(Eigen::VectorXd::Ones(robot.dimv())));
  Eigen::VectorXd violation = Eigen::VectorXd::Zero(robot.dimv());
  EXPECT_TRUE(violation.isApprox(Eigen::VectorXd::Zero(robot.dimv())));
}


TEST_F(FixedBaseRobotTest, numContacts) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  unsigned int num_contact = 5;
  EXPECT_EQ(robot.dimf(), 0);
  for (int i=0; i<num_contact; ++i) {
    robot.addPointContact(2*i+2, 0, 0);
    EXPECT_EQ(robot.dimf(), (i+1)*3);
  }
  for (int i=0; i<num_contact; ++i) {
    robot.removePointContact(2*i+1);
    EXPECT_EQ(robot.dimf(), num_contact*3);
  }
  for (int i=0; i<num_contact; ++i) {
    robot.removePointContact(2*i+2);
    EXPECT_EQ(robot.dimf(), (num_contact-i-1)*3);
  }
}

} // namespace invdynocp 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}