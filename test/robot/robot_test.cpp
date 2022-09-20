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

#include "robotoc/robot/point_contact.hpp"
#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"

#include "urdf_factory.hpp"


namespace robotoc {

class RobotTest : public ::testing::TestWithParam<RobotModelInfo> {
protected:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
  }

  virtual void TearDown() {
  }
};


TEST_P(RobotTest, constructor) {
  Robot robot_empty;
  EXPECT_EQ(robot_empty.dimq(), 0);
  EXPECT_EQ(robot_empty.dimv(), 0);
  EXPECT_EQ(robot_empty.dimu(), 0);
  EXPECT_EQ(robot_empty.max_dimf(), 0);
  EXPECT_EQ(robot_empty.dim_passive(), 0);
  EXPECT_EQ(robot_empty.hasFloatingBase(), false);
  EXPECT_EQ(robot_empty.maxNumContacts(), 0);
  EXPECT_EQ(robot_empty.maxNumPointContacts(), 0);
  EXPECT_EQ(robot_empty.maxNumSurfaceContacts(), 0);
  EXPECT_FALSE(robot_empty.hasFloatingBase());

  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model_ref;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model_ref);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model_ref);
  }
  EXPECT_EQ(robot.dimq(), model_ref.nq);
  EXPECT_EQ(robot.dimv(), model_ref.nv);
  EXPECT_EQ(robot.dimu(), model_ref.nv-robot.dim_passive());
  if (robot.hasFloatingBase()) {
    EXPECT_EQ(robot.dim_passive(), 6);
  }
  else {
    EXPECT_EQ(robot.dim_passive(), 0);
  }
  for (int i=0; i<model_info.point_contacts.size(); ++i) {
    EXPECT_EQ(robot.contactType(i), ContactType::PointContact);
    EXPECT_EQ(robot.contactTypes()[i], ContactType::PointContact);
  }
  for (int i=0; i<model_info.surface_contacts.size(); ++i) {
    EXPECT_EQ(robot.contactType(i+model_info.point_contacts.size()), ContactType::SurfaceContact);
    EXPECT_EQ(robot.contactTypes()[i+model_info.point_contacts.size()], ContactType::SurfaceContact);
  }
  EXPECT_EQ(robot.max_dimf(), 3*model_info.point_contacts.size()+6*model_info.surface_contacts.size());
  EXPECT_EQ(robot.maxNumContacts(), model_info.point_contacts.size()+model_info.surface_contacts.size());
  EXPECT_EQ(robot.maxNumPointContacts(), model_info.point_contacts.size());
  EXPECT_EQ(robot.maxNumSurfaceContacts(), model_info.surface_contacts.size());
  EXPECT_NO_THROW(
    std::cout << robot << std::endl;
  );
  EXPECT_EQ(robot.jointEffortLimit().size(), robot.dimu());
  EXPECT_EQ(robot.jointVelocityLimit().size(), robot.dimu());
  EXPECT_EQ(robot.lowerJointPositionLimit().size(), robot.dimu());
  EXPECT_EQ(robot.upperJointPositionLimit().size(), robot.dimu());
  const Eigen::VectorXd effort_limit = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd velocity_limit = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd lower_position_limit = Eigen::VectorXd::Random(robot.dimu());
  const Eigen::VectorXd upper_position_limit = Eigen::VectorXd::Random(robot.dimu());
  robot.setJointEffortLimit(effort_limit);
  robot.setJointVelocityLimit(velocity_limit);
  robot.setLowerJointPositionLimit(lower_position_limit);
  robot.setUpperJointPositionLimit(upper_position_limit);
  EXPECT_TRUE(robot.jointEffortLimit().isApprox(effort_limit));
  EXPECT_TRUE(robot.jointVelocityLimit().isApprox(velocity_limit));
  EXPECT_TRUE(robot.lowerJointPositionLimit().isApprox(lower_position_limit));
  EXPECT_TRUE(robot.upperJointPositionLimit().isApprox(upper_position_limit));
  pinocchio::Data data = pinocchio::Data(model_ref);
  const double mass_ref = pinocchio::computeTotalMass(model_ref);
  const double weight_ref = - pinocchio::computeTotalMass(model_ref) * pinocchio::Model::gravity981[2];
  EXPECT_DOUBLE_EQ(mass_ref, robot.totalMass());
  EXPECT_DOUBLE_EQ(weight_ref, robot.totalWeight());
  for (int i=0; i<model_ref.frames.size(); ++i) {
    EXPECT_EQ(robot.frameId(robot.frameName(i)), i);
  }
  const auto contact_status = robot.createContactStatus();
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    EXPECT_EQ(robot.contactFrameNames()[i], contact_status.contactFrameNames()[i]);
    EXPECT_EQ(robot.contactFrameNames()[i], contact_status.contactFrameName(i));
  }
}


TEST_P(RobotTest, integrateConfiguration) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  const double integration_length = Eigen::VectorXd::Random(1)[0];
  Eigen::VectorXd q_integrated(robot.dimq()), q_integrated_ref(robot.dimq());
  robot.integrateConfiguration(q, v, integration_length, q_integrated);
  pinocchio::integrate(model, q, integration_length*v, q_integrated_ref);
  EXPECT_TRUE(q_integrated.isApprox(q_integrated_ref));
  Eigen::VectorXd q_integrated2 = q;
  robot.integrateConfiguration(v, integration_length, q_integrated2);
  EXPECT_TRUE(q_integrated2.isApprox(q_integrated));
}


TEST_P(RobotTest, dIntegrateTransport) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const int dim = std::floor(model.nv / 2);
  const Eigen::MatrixXd Jin = Eigen::MatrixXd::Random(dim, model.nv);
  Eigen::MatrixXd dIntdq     = Eigen::MatrixXd::Zero(dim, model.nv);
  Eigen::MatrixXd dIntdv     = Eigen::MatrixXd::Zero(dim, model.nv);
  Eigen::MatrixXd dIntdq_ref = Eigen::MatrixXd::Zero(dim, model.nv);
  Eigen::MatrixXd dIntdv_ref = Eigen::MatrixXd::Zero(dim, model.nv);
  robot.dIntegrateTransport_dq(q, v, Jin, dIntdq);
  robot.dIntegrateTransport_dv(q, v, Jin, dIntdv);
  pinocchio::dIntegrateTransport(model, q, v, Jin.transpose(), dIntdq_ref.transpose(), pinocchio::ARG0);
  pinocchio::dIntegrateTransport(model, q, v, Jin.transpose(), dIntdv_ref.transpose(), pinocchio::ARG1);
  EXPECT_TRUE(dIntdq.isApprox(dIntdq_ref));
  EXPECT_TRUE(dIntdv.isApprox(dIntdv_ref));
}


TEST_P(RobotTest, subtractConfiguration) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }

  const Eigen::VectorXd q_plus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd q_minus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));    
  Eigen::VectorXd q_diff(model.nv), q_diff_ref(model.nv);
  robot.subtractConfiguration(q_plus, q_minus, q_diff);
  pinocchio::difference(model, q_minus, q_plus, q_diff_ref);
  EXPECT_TRUE(q_diff.isApprox(q_diff_ref));
}


TEST_P(RobotTest, subtractConfigurationDerivatives) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }

  const Eigen::VectorXd q_plus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd q_minus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  Eigen::MatrixXd dSubdq_plus =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_plus_ref =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus_ref =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  robot.dSubtractConfiguration_dqf(q_plus, q_minus, dSubdq_plus);
  robot.dSubtractConfiguration_dq0(q_plus, q_minus, dSubdq_minus);
  pinocchio::dDifference(model, q_minus, q_plus, dSubdq_plus_ref, pinocchio::ARG1);
  pinocchio::dDifference(model, q_minus, q_plus, dSubdq_minus_ref, pinocchio::ARG0);
  EXPECT_TRUE(dSubdq_plus.isApprox(dSubdq_plus_ref));
  EXPECT_TRUE(dSubdq_minus.isApprox(dSubdq_minus_ref));
  Eigen::MatrixXd dSubdq_plus_inv =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus_inv =  Eigen::MatrixXd::Zero(model.nv, model.nv);
}


TEST_P(RobotTest, frameKinematics) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  if (model_info.point_contacts.empty() || model_info.surface_contacts.empty()) {
    return;
  }

  int frame_id = 0;
  if (model_info.point_contacts.empty()) {
    frame_id = robot.frameId(model_info.surface_contacts[0].frame);
  }
  else {
    frame_id = robot.frameId(model_info.point_contacts[0].frame);
  }

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  robot.updateKinematics(q, v);
  const auto frame_position = robot.framePosition(frame_id);
  const auto frame_rotation = robot.frameRotation(frame_id);
  const auto frame_placement = robot.framePlacement(frame_id);
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, model.nv);
  robot.getFrameJacobian(frame_id, J);
  pinocchio::forwardKinematics(model, data, q);
  pinocchio::updateFramePlacements(model, data);
  const auto frame_position_ref = data.oMf[frame_id].translation();
  EXPECT_TRUE(frame_position.isApprox(frame_position_ref));
  const auto frame_rotation_ref = data.oMf[frame_id].rotation();
  EXPECT_TRUE(frame_rotation.isApprox(frame_rotation_ref));
  const auto frame_placement_ref = data.oMf[frame_id];
  EXPECT_TRUE(frame_placement.isApprox(frame_placement_ref));
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeJointJacobians(model, data, q);
  Eigen::MatrixXd J_ref = Eigen::MatrixXd::Zero(6, model.nv);
  pinocchio::getFrameJacobian(model, data, frame_id, pinocchio::LOCAL, J_ref);
  EXPECT_TRUE(J.isApprox(J_ref));
  pinocchio::centerOfMass(model, data, q, true);
  EXPECT_TRUE(data.com[0].isApprox(robot.CoM()));
  pinocchio::jacobianCenterOfMass(model, data, true);
  Eigen::MatrixXd Jcom = Eigen::MatrixXd::Zero(3, model.nv);
  robot.getCoMJacobian(Jcom);
  EXPECT_TRUE(Jcom.isApprox(data.Jcom));
}


TEST_P(RobotTest, transformFromLocalToWorld) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }

  if (model_info.point_contacts.empty() || model_info.surface_contacts.empty()) {
    return;
  }

  int frame_id = 0;
  if (model_info.point_contacts.empty()) {
    frame_id = robot.frameId(model_info.surface_contacts[0].frame);
  }
  else {
    frame_id = robot.frameId(model_info.point_contacts[0].frame);
  }

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  robot.updateKinematics(q, v);
  const auto frame_rotation = robot.frameRotation(frame_id);
  const Eigen::Vector3d vec_local = Eigen::Vector3d::Random();
  const Eigen::Vector3d vec_world_ref = frame_rotation * vec_local;
  Eigen::Vector3d vec_world = Eigen::Vector3d::Zero();
  robot.transformFromLocalToWorld(frame_id, vec_local, vec_world);
  EXPECT_TRUE(vec_world.isApprox(vec_world_ref));
  Eigen::MatrixXd J_ref = Eigen::MatrixXd::Zero(6, model.nv);
  robot.getFrameJacobian(frame_id, J_ref);
  for (int i=0; i<model.nv; ++i) {
    J_ref.template topRows<3>().col(i) 
        = J_ref.template bottomRows<3>().col(i).cross(vec_world);
  }
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, model.nv);
  robot.getJacobianTransformFromLocalToWorld(frame_id, vec_world, J);
  EXPECT_TRUE(J.isApprox(J_ref));
}


TEST_P(RobotTest, baumgarte) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  std::vector<PointContact> point_contacts_ref; 
  std::vector<SurfaceContact> surface_contacts_ref; 
  for (const auto& e : model_info.point_contacts) {
    point_contacts_ref.push_back(PointContact(model, e));
  }
  for (const auto& e : model_info.surface_contacts) {
    surface_contacts_ref.push_back(SurfaceContact(model, e));
  }

  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  robot.updateKinematics(q, v, a);
  robot.computeBaumgarteResidual(contact_status, residual.head(contact_status.dimf()));
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  int dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_ref[i].computeBaumgarteResidual(
          model, data, contact_status.contactPosition(i), 
          residual_ref.segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeBaumgarteResidual(
          model, data, contact_status.contactPlacement(i+point_contacts_ref.size()), 
          residual_ref.segment<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));

  Eigen::MatrixXd baumgarte_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeBaumgarteDerivatives(contact_status, 
                                    baumgarte_partial_q.topRows(contact_status.dimf()), 
                                    baumgarte_partial_v.topRows(contact_status.dimf()), 
                                    baumgarte_partial_a.topRows(contact_status.dimf()));
  Eigen::MatrixXd baumgarte_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      point_contacts_ref[i].computeBaumgarteDerivatives(
          model, data, baumgarte_partial_q_ref.middleRows<3>(dimf), 
          baumgarte_partial_v_ref.middleRows<3>(dimf), 
          baumgarte_partial_a_ref.middleRows<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeBaumgarteDerivatives(
          model, data, baumgarte_partial_q_ref.middleRows<6>(dimf), 
          baumgarte_partial_v_ref.middleRows<6>(dimf), 
          baumgarte_partial_a_ref.middleRows<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(baumgarte_partial_q.isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(baumgarte_partial_v.isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(baumgarte_partial_a.isApprox(baumgarte_partial_a_ref));
}


TEST_P(RobotTest, impulseVelocity) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  std::vector<PointContact> point_contacts_ref; 
  std::vector<SurfaceContact> surface_contacts_ref; 
  for (const auto& e : model_info.point_contacts) {
    point_contacts_ref.push_back(PointContact(model, e));
  }
  for (const auto& e : model_info.surface_contacts) {
    surface_contacts_ref.push_back(SurfaceContact(model, e));
  }

  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  robot.updateKinematics(q, v);
  robot.computeImpulseVelocityResidual(impulse_status, residual.head(impulse_status.dimi()));
  pinocchio::forwardKinematics(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  int dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_ref[i].computeContactVelocityResidual(
          model, data, residual_ref.segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeContactVelocityResidual(
          model, data, residual_ref.segment<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));

  Eigen::MatrixXd velocity_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd velocity_partial_v = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeImpulseVelocityDerivatives(impulse_status, velocity_partial_q.topRows(impulse_status.dimi()), 
                                          velocity_partial_v.topRows(impulse_status.dimi()));
  Eigen::MatrixXd velocity_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd velocity_partial_v_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_ref[i].computeContactVelocityDerivatives(
          model, data, 
          velocity_partial_q_ref.middleRows<3>(dimf), 
          velocity_partial_v_ref.middleRows<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeContactVelocityDerivatives(
          model, data, 
          velocity_partial_q_ref.middleRows<6>(dimf), 
          velocity_partial_v_ref.middleRows<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(velocity_partial_q.isApprox(velocity_partial_q_ref));
  EXPECT_TRUE(velocity_partial_v.isApprox(velocity_partial_v_ref));
}


TEST_P(RobotTest, contactPosition) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  std::vector<PointContact> point_contacts_ref; 
  std::vector<SurfaceContact> surface_contacts_ref; 
  for (const auto& e : model_info.point_contacts) {
    point_contacts_ref.push_back(PointContact(model, e));
  }
  for (const auto& e : model_info.surface_contacts) {
    surface_contacts_ref.push_back(SurfaceContact(model, e));
  }

  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  robot.updateKinematics(q, v);
  robot.computeContactPositionResidual(impulse_status, residual.head(impulse_status.dimi()));
  pinocchio::forwardKinematics(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  int dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_ref[i].computeContactPositionResidual(
          model, data, impulse_status.contactPosition(i), 
          residual_ref.segment<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeContactPositionResidual(
          model, data, impulse_status.contactPlacement(i+point_contacts_ref.size()), 
          residual_ref.segment<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));

  Eigen::MatrixXd contact_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeContactPositionDerivative(impulse_status, contact_partial_q.topRows(impulse_status.dimi()));

  Eigen::MatrixXd contact_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  dimf = 0;
  for (int i=0; i<point_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      point_contacts_ref[i].computeContactPositionDerivative(
          model, data, contact_partial_q_ref.middleRows<3>(dimf));
      dimf += 3;
    }
  }
  for (int i=0; i<surface_contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i+point_contacts_ref.size())) {
      surface_contacts_ref[i].computeContactPositionDerivative(
          model, data, contact_partial_q_ref.middleRows<6>(dimf));
      dimf += 6;
    }
  }
  EXPECT_TRUE(contact_partial_q.isApprox(contact_partial_q_ref));
}


TEST_P(RobotTest, RNEA) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  std::vector<Vector6d> f;
  for (const auto& e : model_info.point_contacts) {
    f.push_back(Vector6d::Random());
  }
  for (const auto& e : model_info.surface_contacts) {
    f.push_back(Vector6d::Random());
  }
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  robot.setContactForces(contact_status, f);

  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
  robot.RNEA(q, v, a, tau);

  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(model.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<model_info.point_contacts.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      PointContact(model, model_info.point_contacts[i]).computeJointForceFromContactForce(
          f[i].template head<3>(), fjoint);
    }
  }
  for (int i=0; i<model_info.surface_contacts.size(); ++i) {
    if (contact_status.isContactActive(i+model_info.point_contacts.size())) {
      SurfaceContact(model, model_info.surface_contacts[i]).computeJointForceFromContactWrench(
          f[i+model_info.point_contacts.size()], fjoint);
    }
  }
  Eigen::VectorXd tau_ref = pinocchio::rnea(model, data, q, v, a, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));

  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = dRNEA_da;
  robot.RNEADerivatives(q, v, a, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::computeRNEADerivatives(model, data, q, v, a, fjoint, dRNEA_dq_ref, dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
}


TEST_P(RobotTest, RNEAImpulse) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  model.gravity.linear().setZero();
  auto data = pinocchio::Data(model);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(model.nv);
  std::vector<Vector6d> f;
  for (const auto& e : model_info.point_contacts) {
    f.push_back(Vector6d::Random());
  }
  for (const auto& e : model_info.surface_contacts) {
    f.push_back(Vector6d::Random());
  }
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  robot.setImpulseForces(impulse_status, f);

  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
  robot.RNEAImpulse(q, dv, tau);

  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(model.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<model_info.point_contacts.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      PointContact(model, model_info.point_contacts[i]).computeJointForceFromContactForce(
          f[i].template head<3>(), fjoint);
    }
  }
  for (int i=0; i<model_info.surface_contacts.size(); ++i) {
    if (impulse_status.isImpulseActive(i+model_info.point_contacts.size())) {
      SurfaceContact(model, model_info.surface_contacts[i]).computeJointForceFromContactWrench(
          f[i+model_info.point_contacts.size()], fjoint);
    }
  }
  Eigen::VectorXd tau_ref = pinocchio::rnea(model, data, q, Eigen::VectorXd::Zero(model.nv), dv, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));

  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_ddv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_ddv_ref = dRNEA_ddv;
  Eigen::MatrixXd dRNEA_dv_ref = Eigen::MatrixXd::Zero(model.nv, model.nv);
  robot.RNEAImpulseDerivatives(q, dv, dRNEA_dq, dRNEA_ddv);
  pinocchio::computeRNEADerivatives(model, data, q, Eigen::VectorXd::Zero(model.nv), dv, fjoint, 
                                    dRNEA_dq_ref, dRNEA_dv_ref, dRNEA_ddv_ref);
  dRNEA_ddv_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_ddv_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_ddv.isApprox(dRNEA_ddv_ref));
}


TEST_P(RobotTest, MJtJinv) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  pinocchio::Model model;
  if (model_info.base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(model_info.urdf_path, pinocchio::JointModelFreeFlyer(), model);
  }
  else {
    pinocchio::urdf::buildModel(model_info.urdf_path, model);
  }
  auto data = pinocchio::Data(model);

  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  const int dimf = contact_status.dimf();
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(model.nv, model.nv);
  robot.RNEADerivatives(q, v, a, dRNEA_dq, dRNEA_dv, dRNEA_da);
  if (dimf > 0) {
    Eigen::MatrixXd baumgarte_partial_q = Eigen::MatrixXd::Zero(dimf, model.nv);
    Eigen::MatrixXd baumgarte_partial_v = Eigen::MatrixXd::Zero(dimf, model.nv);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(dimf, model.nv);
    robot.computeBaumgarteDerivatives(contact_status, baumgarte_partial_q, 
                                      baumgarte_partial_v, J);
    Eigen::MatrixXd MJtJinv = Eigen::MatrixXd::Zero(model.nv+dimf, model.nv+dimf);
    robot.computeMJtJinv(dRNEA_da, J, MJtJinv);
    Eigen::MatrixXd MJtJ = Eigen::MatrixXd::Zero(model.nv+dimf, model.nv+dimf);
    MJtJ.topLeftCorner(model.nv, model.nv) = dRNEA_da;
    MJtJ.topRightCorner(model.nv, dimf) = J.transpose();
    MJtJ.bottomLeftCorner(dimf, model.nv) = J;
    MJtJ.bottomRightCorner(dimf, dimf).diagonal().array() = - model_info.contact_inv_damping;
    const Eigen::MatrixXd MJtJinv_ref = MJtJ.inverse();
    EXPECT_TRUE(MJtJinv.isApprox(MJtJinv_ref));
  }
  Eigen::MatrixXd Minv = dRNEA_da;
  robot.computeMinv(dRNEA_da, Minv);
  EXPECT_TRUE((dRNEA_da*Minv).isIdentity());
  EXPECT_TRUE((Minv*dRNEA_da).isIdentity());
}


TEST_P(RobotTest, generateFeasibleConfiguration) {
  const auto model_info = GetParam();
  Robot robot(model_info);
  const Eigen::VectorXd q = robot.generateFeasibleConfiguration();
  EXPECT_TRUE(q.size() == robot.dimq());
  Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  if (robot.hasFloatingBase()) {
    for (int i=0; i<robot.dimq()-robot.dim_passive()-1; ++i) {
      EXPECT_TRUE(q(robot.dim_passive()+1+i) >= qmin(i));
      EXPECT_TRUE(q(robot.dim_passive()+1+i) <= qmax(i));
    }
  }
  else {
    for (int i=0; i<robot.dimq(); ++i) {
      EXPECT_TRUE(q(i) >= qmin(i));
      EXPECT_TRUE(q(i) <= qmax(i));
    }
  }
  Eigen::VectorXd q_norm  = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_norm);
  EXPECT_TRUE(q_norm.size() == robot.dimq());
}


auto manipulatorInfo = [](const bool contact) {
  auto info = RobotModelInfo::Manipulator(testhelper::RobotManipulatorURDF());
  if (contact) {
    info.point_contacts.push_back(ContactModelInfo("iiwa_link_ee_kuka", 0.1));
  }
  return info;
};


auto quadrupedInfo = [](const bool contact) {
  std::vector<ContactModelInfo> point_contacts;
  if (contact) {
    point_contacts.push_back(ContactModelInfo("LF_FOOT", 0.1));
    point_contacts.push_back(ContactModelInfo("LH_FOOT", 0.1));
    point_contacts.push_back(ContactModelInfo("RF_FOOT", 0.1));
    point_contacts.push_back(ContactModelInfo("RH_FOOT", 0.1));
  }
  return RobotModelInfo::Quadruped(testhelper::QuadrupedURDF(), point_contacts);
};

auto humanoidInfo = [](const bool contact) {
  std::vector<ContactModelInfo> surface_contacts;
  if (contact) {
    surface_contacts.push_back(ContactModelInfo("l_sole", 0.1));
    surface_contacts.push_back(ContactModelInfo("r_sole", 0.1));
  }
  return RobotModelInfo::Humanoid(testhelper::HumanoidURDF(), surface_contacts);
};


INSTANTIATE_TEST_SUITE_P(
  TestWithMultipleRobots, RobotTest, 
  ::testing::Values(manipulatorInfo(false), 
                    manipulatorInfo(true),
                    quadrupedInfo(false),
                    quadrupedInfo(true),
                    humanoidInfo(false),
                    humanoidInfo(true))
);

// INSTANTIATE_TEST_SUITE_P(
//   TestWithMultipleRobots, RobotTest, 
//   ::testing::Values(manipulatorInfo(false), 
//                     manipulatorInfo(true),
//                     quadrupedInfo(false),
//                     quadrupedInfo(true))
// );

} // namespace robotoc 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}