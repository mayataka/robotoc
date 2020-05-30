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

#include "common/memory_manager.hpp"
#include "robot/point_contact.hpp"
#include "robot/robot.hpp"


namespace robotcgmres {

class FixedBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    contact_frame_id_ = 18;
    max_point_contacts_ = 1;
    q_ = memorymanager::NewVector(dimq_);
    v_ = memorymanager::NewVector(dimq_);
    a_ = memorymanager::NewVector(dimq_);
    Eigen::Map<Eigen::VectorXd>(q_, dimq_) 
        = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                         Eigen::VectorXd::Ones(dimq_));
    Eigen::Map<Eigen::VectorXd>(v_, dimq_) = Eigen::VectorXd::Random(dimq_);
    Eigen::Map<Eigen::VectorXd>(a_, dimq_) = Eigen::VectorXd::Random(dimq_);
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[0];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[0];
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(a_);
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, contact_frame_id_, max_point_contacts_;
  double *q_, *v_, *a_;
  double baumgarte_alpha_, baumgarte_beta_;
};


TEST_F(FixedBaseRobotTest, dim) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimq_);
  EXPECT_EQ(robot.dimfmax(), 3*max_point_contacts_);
  EXPECT_EQ(robot.dimf(), 0);
  EXPECT_EQ(robot.dim_passive(), 0);
}


TEST_F(FixedBaseRobotTest, integrateConfiguration) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *q_plus = memorymanager::NewVector(dimq_);
  const double integration_length = std::abs(Eigen::VectorXd::Random(2)[1]);
  robot.integrateConfiguration(q_, v_, integration_length, q_plus);
  double *q_plus_ref = memorymanager::NewVector(dimq_);
  for (int i=0; i<dimq_; ++i) {
    q_plus_ref[i] = q_[i] + integration_length * v_[i];
    EXPECT_DOUBLE_EQ(q_plus[i], q_plus_ref[i]);
  }
  robot.integrateConfiguration(v_, integration_length, q_);
  for (int i=0; i<dimq_; ++i) {
    q_plus_ref[i] = q_[i];
  }
  memorymanager::DeleteVector(q_plus);
  memorymanager::DeleteVector(q_plus_ref);
}


TEST_F(FixedBaseRobotTest, baumgarteResidualAndDerivatives) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  double *residual = memorymanager::NewVector(3*1);
  robot.updateKinematics(q_, v_, a_);
  robot.computeBaumgarteResidual(residual);
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q_, dimq_),
                               Eigen::Map<const Eigen::VectorXd>(v_, dimq_),
                               Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  pinocchio::computeForwardKinematicsDerivatives(
      model_, data_, Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  pinocchio::updateFramePlacements(model_, data_);
  contact_ref.resetContactPointByCurrentKinematics(data_);
  double *residual_ref = memorymanager::NewVector(3*1);
  contact_ref.computeBaumgarteResidual(model_, data_, residual_ref);
  EXPECT_TRUE(Eigen::Map<Eigen::Vector3d>(residual).isApprox(Eigen::Map<Eigen::Vector3d>(residual_ref)));
  memorymanager::DeleteVector(residual);
  memorymanager::DeleteVector(residual_ref);
  double *dBaum_dq_dot_vec = memorymanager::NewVector(dimq_);
  double *dBaum_dv_dot_vec = memorymanager::NewVector(dimq_);
  double *dBaum_da_dot_vec = memorymanager::NewVector(dimq_);
  double *vec = memorymanager::NewVector(3*1);
  Eigen::Map<Eigen::VectorXd>(vec, 3) = Eigen::VectorXd::Random(3*1);
  robot.computeBaumgarteDerivativesDotVec(vec, dBaum_dq_dot_vec, 
                                          dBaum_dv_dot_vec, dBaum_da_dot_vec);
  Eigen::MatrixXd dBaum_dq, dBaum_dv, dBaum_da;
  dBaum_dq = Eigen::MatrixXd::Zero(3*1, dimq_);
  dBaum_dv = Eigen::MatrixXd::Zero(3*1, dimq_);
  dBaum_da = Eigen::MatrixXd::Zero(3*1, dimq_);
  contact_ref.computeBaumgarteDerivatives(model_, data_, dBaum_dq, dBaum_dv, 
                                          dBaum_da, 0, 0);
  Eigen::VectorXd dBaum_dq_dot_vec_ref
      = dBaum_dq.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, 3*1);
  Eigen::VectorXd dBaum_dv_dot_vec_ref
      = dBaum_dv.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, 3*1);
  Eigen::VectorXd dBaum_da_dot_vec_ref
      = dBaum_da.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, 3*1);
  EXPECT_TRUE(dBaum_dq_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dBaum_dq_dot_vec, dimq_)));
  EXPECT_TRUE(dBaum_dv_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dBaum_dv_dot_vec, dimq_)));
  EXPECT_TRUE(dBaum_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dBaum_da_dot_vec, dimq_)));
  memorymanager::DeleteVector(dBaum_dq_dot_vec);
  memorymanager::DeleteVector(dBaum_dv_dot_vec);
  memorymanager::DeleteVector(dBaum_da_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, RNEAWithoutFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *tau = memorymanager::NewVector(dimq_);
  robot.RNEA(q_, v_, a_, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  EXPECT_TRUE(tau_ref.isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  memorymanager::DeleteVector(tau);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesTransDotVecWithoutFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *dRNEA_dq_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec = memorymanager::NewVector(dimq_);
  robot.RNEADerivativesTransDotVec(q_, v_, a_, vec, dRNEA_dq_dot_vec, 
                                   dRNEA_dv_dot_vec, dRNEA_da_dot_vec);
  pinocchio::computeRNEADerivatives(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd dRNEA_dq_dot_vec_ref 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_dv_dot_vec_ref 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_da_dot_vec_ref 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_dq_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dv_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_)));
  memorymanager::DeleteVector(dRNEA_dq_dot_vec);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec);
  memorymanager::DeleteVector(dRNEA_da_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesWithAccTransDotVecWithoutFext) {
  Robot robot(urdf_file_name_, 0);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *result = memorymanager::NewVector(dimq_);
  robot.RNEADerivativesTransDotVec(q_, vec, result);
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd result_ref
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(result_ref.isApprox(Eigen::Map<Eigen::VectorXd>(result, dimq_)));
  memorymanager::DeleteVector(vec);
  memorymanager::DeleteVector(result);
}


TEST_F(FixedBaseRobotTest, RNEAWithFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *tau = memorymanager::NewVector(dimq_);
  double *tau0 = memorymanager::NewVector(dimq_);
  double *fext = memorymanager::NewVector(3*1);
  Eigen::Map<Eigen::VectorXd>(fext, 3*1) = Eigen::VectorXd::Random(3*1);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  EXPECT_EQ(robot.dimf(), 3*1);
  robot.RNEA(q_, v_, a_, tau0);
  robot.RNEA(q_, v_, a_, fext, tau);
  EXPECT_TRUE((Eigen::Map<Eigen::VectorXd>(tau0, dimq_)).isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  robot.setFext(fext);
  robot.RNEA(q_, v_, a_, fext, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  Eigen::VectorXd tau_ref = pinocchio::rnea(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_),
      fjoint);
  EXPECT_TRUE(tau_ref.isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  memorymanager::DeleteVector(tau);
  memorymanager::DeleteVector(tau0);
  memorymanager::DeleteVector(fext);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesTransDotVecWithFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *dRNEA_dq_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dfext_dot_vec = memorymanager::NewVector(3*1);
  double *dRNEA_dq_dot_vec0 = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec0 = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec0 = memorymanager::NewVector(dimq_);
  double *fext = memorymanager::NewVector(3*1);
  Eigen::Map<Eigen::VectorXd>(fext, 3*1) = Eigen::VectorXd::Random(3*1);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  robot.RNEADerivativesTransDotVec(q_, v_, a_,  vec, dRNEA_dq_dot_vec0, 
                                   dRNEA_dv_dot_vec0, dRNEA_da_dot_vec0);
  robot.RNEADerivativesTransDotVec(q_, v_, a_, fext, vec, dRNEA_dq_dot_vec, 
                                   dRNEA_dv_dot_vec, dRNEA_da_dot_vec);
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec0, dimq_)));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec0, dimq_)));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec0, dimq_)));
  robot.setFext(fext);
  robot.RNEADerivativesTransDotVec(q_, v_, a_, fext, vec, dRNEA_dq_dot_vec, 
                                   dRNEA_dv_dot_vec, dRNEA_da_dot_vec);
  robot.updateKinematics(q_, v_, a_);
  robot.dRNEAdfextTransDotVec(vec, dRNEA_dfext_dot_vec);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.computeJointForceFromContactForce(fext, fjoint);
  pinocchio::computeRNEADerivatives(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_), fjoint);
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd dRNEA_dq_dot_vec_ref 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_dv_dot_vec_ref 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_da_dot_vec_ref 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::MatrixXd dtau_dfext = Eigen::MatrixXd::Zero(6, dimq_);
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::computeJointJacobians(model_, data_, 
                                   Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::updateFramePlacements(model_, data_);
  contact_ref.computeContactJacobian(model_, data_, dtau_dfext, 0, 0);
  Eigen::VectorXd dRNEA_dfext_dot_vec_ref 
      = dtau_dfext.topRows(3) * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_dq_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dv_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dfext_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dfext_dot_vec, 3)));
  memorymanager::DeleteVector(fext);
  memorymanager::DeleteVector(dRNEA_dq_dot_vec0);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec0);
  memorymanager::DeleteVector(dRNEA_da_dot_vec0);
  memorymanager::DeleteVector(dRNEA_dq_dot_vec);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec);
  memorymanager::DeleteVector(dRNEA_da_dot_vec);
  memorymanager::DeleteVector(dRNEA_dfext_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesWithAccAndFextTransDotVecWithFext) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  robot.updateKinematics(q_, v_, a_);
  double *dRNEA_da_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dfext_dot_vec = memorymanager::NewVector(3*1);
  robot.RNEADerivativesTransDotVec(q_, vec, dRNEA_da_dot_vec);
  robot.dRNEAdfextTransDotVec(vec, dRNEA_dfext_dot_vec);
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::computeJointJacobians(model_, data_, 
                                   Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::updateFramePlacements(model_, data_);
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd dRNEA_da_dot_vec_ref
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_)));
  PointContact contact_ref(model_, contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  Eigen::MatrixXd dtau_dfext = Eigen::MatrixXd::Zero(6, dimq_);
  contact_ref.computeContactJacobian(model_, data_, dtau_dfext, 0, 0);
  Eigen::VectorXd dRNEA_dfext_dot_vec_ref 
      = dtau_dfext.topRows(3) * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_dfext_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dfext_dot_vec, 3*1)));
  memorymanager::DeleteVector(dRNEA_da_dot_vec);
  memorymanager::DeleteVector(dRNEA_dfext_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, passiveJoints) {
  Robot robot(urdf_file_name_, 0);
  double* tau = memorymanager::NewVector(robot.dimv());
  for (int i=0; i<robot.dimv(); ++i) {
    tau[i] = 1.0;
  }
  robot.setPassiveTorques(tau);
  for (int i=0; i<robot.dimv(); ++i) {
    EXPECT_DOUBLE_EQ(tau[i], 1.0);
  }
  double* residual = memorymanager::NewVector(robot.dimv());
  robot.passiveConstraintsResidual(tau, residual);
  for (int i=0; i<robot.dimv(); ++i) {
    EXPECT_DOUBLE_EQ(residual[i], 0.0);
  }
  robot.addPassiveConstraintsDerivativeDotVec(tau, residual);
  for (int i=0; i<robot.dimv(); ++i) {
    EXPECT_DOUBLE_EQ(residual[i], 0.0);
  }
  memorymanager::DeleteVector(residual);
  memorymanager::DeleteVector(tau);
}


TEST_F(FixedBaseRobotTest, numContacts) {
  Robot robot(urdf_file_name_, max_point_contacts_);
  int num_contact = 5;
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

} // namespace robot_cgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}