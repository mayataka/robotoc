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
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

class RobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    pinocchio::urdf::buildModel(fixed_base_urdf, fixed_base_robot);
    pinocchio::urdf::buildModel(floating_base_urdf, 
                                pinocchio::JointModelFreeFlyer(), 
                                floating_base_robot);
    fixed_base_data = pinocchio::Data(fixed_base_robot);
    floating_base_data = pinocchio::Data(floating_base_robot);
    fixed_base_contact_frames = {18};
    floating_base_contact_frames = {12, 22, 32, 42};
    baumgarte_weight_on_velocity = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
    baumgarte_weight_on_position = 10 * std::abs(Eigen::VectorXd::Random(1)[0]);
    baumgarte_weights.first = baumgarte_weight_on_velocity;
    baumgarte_weights.second = baumgarte_weight_on_position;
  }

  virtual void TearDown() {
  }

  void testConstructorAndSetter(const std::string& path_to_urdf, 
                                const BaseJointType& base_joint_type,
                                const std::vector<int>& contact_frames) const;
  void testIntegrateConfiguration(const std::string& path_to_urdf, 
                                  const BaseJointType& base_joint_type,
                                  pinocchio::Model& model) const;
  void testIntegrateConfigurationDerivatives(const std::string& path_to_urdf, 
                                             const BaseJointType& base_joint_type,
                                             pinocchio::Model& model) const;
  void testSubtractConfiguration(const std::string& path_to_urdf, 
                                 const BaseJointType& base_joint_type,
                                 pinocchio::Model& model) const;
  void testSubtractConfigurationDerivatives(const std::string& path_to_urdf, 
                                            const BaseJointType& base_joint_type,
                                            pinocchio::Model& model) const;
  void testFrameKinematics(const std::string& path_to_urdf, 
                           const BaseJointType& base_joint_type,
                           pinocchio::Model& model, pinocchio::Data& data, 
                           const int frame_id) const;
  void testBaumgarte(const std::string& path_to_urdf, 
                     const BaseJointType& base_joint_type,
                     pinocchio::Model& model, pinocchio::Data& data, 
                     const std::vector<int>& frames) const;
  void testBaumgarteTimeStep(const std::string& path_to_urdf, 
                             const BaseJointType& base_joint_type,
                             pinocchio::Model& model, pinocchio::Data& data, 
                             const std::vector<int>& frames) const;
  void testImpulseVelocity(const std::string& path_to_urdf, 
                           const BaseJointType& base_joint_type,
                           pinocchio::Model& model, pinocchio::Data& data, 
                           const std::vector<int>& frames) const;
  void testContactPosition(const std::string& path_to_urdf, 
                           const BaseJointType& base_joint_type,
                           pinocchio::Model& model, pinocchio::Data& data, 
                           const std::vector<int>& frames) const;
  void testRNEA(const std::string& path_to_urdf, 
                const BaseJointType& base_joint_type,
                pinocchio::Model& model, pinocchio::Data& data) const;
  void testRNEA(const std::string& path_to_urdf, 
                const BaseJointType& base_joint_type,
                pinocchio::Model& model, pinocchio::Data& data, 
                const std::vector<int>& frames) const;
  void testRNEAImpulse(const std::string& path_to_urdf, 
                       const BaseJointType& base_joint_type, 
                       pinocchio::Model& model, pinocchio::Data& data, 
                       const std::vector<int>& frames) const;
  void testMJtJinv(const std::string& path_to_urdf, 
                   const BaseJointType& base_joint_type, 
                   pinocchio::Model& model, pinocchio::Data& data, 
                   const std::vector<int>& frames) const;
  void testGenConfiguration(const std::string& path_to_urdf, 
                            const BaseJointType& base_joint_type, 
                            pinocchio::Model& model, pinocchio::Data& data) const;

  std::string fixed_base_urdf, floating_base_urdf;
  pinocchio::Model fixed_base_robot, floating_base_robot;
  pinocchio::Data fixed_base_data, floating_base_data;
  std::vector<int> fixed_base_contact_frames, floating_base_contact_frames;
  double baumgarte_weight_on_velocity, baumgarte_weight_on_position;
  std::pair<double, double> baumgarte_weights;
};


void RobotTest::testConstructorAndSetter(const std::string& path_to_urdf, 
                                         const BaseJointType& base_joint_type,
                                         const std::vector<int>& contact_frames) const {
  Robot robot_empty;
  EXPECT_EQ(robot_empty.dimq(), 0);
  EXPECT_EQ(robot_empty.dimv(), 0);
  EXPECT_EQ(robot_empty.dimu(), 0);
  EXPECT_EQ(robot_empty.max_dimf(), 0);
  EXPECT_EQ(robot_empty.dim_passive(), 0);
  EXPECT_EQ(robot_empty.maxPointContacts(), 0);
  EXPECT_FALSE(robot_empty.hasFloatingBase());
  Robot robot(path_to_urdf, base_joint_type);
  pinocchio::Model robot_ref;
  if (base_joint_type == BaseJointType::FloatingBase) {
    pinocchio::urdf::buildModel(path_to_urdf, pinocchio::JointModelFreeFlyer(), 
                                robot_ref);
  }
  else {
    pinocchio::urdf::buildModel(path_to_urdf, robot_ref);
  }
  EXPECT_EQ(robot.dimq(), robot_ref.nq);
  EXPECT_EQ(robot.dimv(), robot_ref.nv);
  EXPECT_EQ(robot.dimu(), robot_ref.nv-robot.dim_passive());
  EXPECT_EQ(robot.max_dimf(), 0);
  if (robot.hasFloatingBase()) {
    EXPECT_EQ(robot.dim_passive(), 6);
  }
  else {
    EXPECT_EQ(robot.dim_passive(), 0);
  }
  EXPECT_EQ(robot.maxPointContacts(), 0);
  Robot robot_contact(path_to_urdf, base_joint_type, contact_frames, baumgarte_weights);
  EXPECT_EQ(robot_contact.dimq(), robot_ref.nq);
  EXPECT_EQ(robot_contact.dimv(), robot_ref.nv);
  EXPECT_EQ(robot_contact.dimu(), robot_ref.nv-robot_contact.dim_passive());
  EXPECT_EQ(robot_contact.max_dimf(), 3*contact_frames.size());
  if (robot.hasFloatingBase()) {
    EXPECT_EQ(robot_contact.dim_passive(), 6);
  }
  else {
    EXPECT_EQ(robot_contact.dim_passive(), 0);
  }
  EXPECT_EQ(robot_contact.maxPointContacts(), contact_frames.size());
  robot_contact.printRobotModel();
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
  pinocchio::Data data = pinocchio::Data(robot_ref);
  const double weight_ref = - pinocchio::computeTotalMass(robot_ref) * pinocchio::Model::gravity981[2];
  EXPECT_DOUBLE_EQ(weight_ref, robot.totalWeight());
}


void RobotTest::testIntegrateConfiguration(const std::string& path_to_urdf, 
                                           const BaseJointType& base_joint_type,
                                           pinocchio::Model& model) const {
  Robot robot(path_to_urdf, base_joint_type);
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


void RobotTest::testIntegrateConfigurationDerivatives(const std::string& path_to_urdf, 
                                                      const BaseJointType& base_joint_type,
                                                      pinocchio::Model& model) const {
  Robot robot(path_to_urdf, base_joint_type);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  Eigen::MatrixXd dIntdq     = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dIntdv     = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dIntdq_ref = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dIntdv_ref = Eigen::MatrixXd::Zero(model.nv, model.nv);
  robot.dIntegratedConfiguration(q, v, dIntdq);
  robot.dIntegratedVelocity(q, v, dIntdv);
  pinocchio::dIntegrate(model, q, v, dIntdq_ref, pinocchio::ARG0);
  pinocchio::dIntegrate(model, q, v, dIntdv_ref, pinocchio::ARG1);
  EXPECT_TRUE(dIntdq.isApprox(dIntdq_ref));
  EXPECT_TRUE(dIntdv.isApprox(dIntdv_ref));
  std::cout << "dIntdq" << std::endl;
  std::cout << dIntdq << std::endl;
  std::cout << "dIntdv" << std::endl;
  std::cout << dIntdv << std::endl;
}


void RobotTest::testSubtractConfiguration(const std::string& path_to_urdf, 
                                          const BaseJointType& base_joint_type,
                                          pinocchio::Model& model) const {
  Robot robot(path_to_urdf, base_joint_type);
  const Eigen::VectorXd q_plus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd q_minus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));    
  Eigen::VectorXd q_diff(model.nv), q_diff_ref(model.nv);
  robot.subtractConfiguration(q_plus, q_minus, q_diff);
  pinocchio::difference(model, q_minus, q_plus, q_diff_ref);
  EXPECT_TRUE(q_diff.isApprox(q_diff_ref));
}


void RobotTest::testSubtractConfigurationDerivatives(const std::string& path_to_urdf, 
                                                     const BaseJointType& base_joint_type,
                                                     pinocchio::Model& model) const {
  Robot robot(path_to_urdf, base_joint_type);
  const Eigen::VectorXd q_plus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd q_minus = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  Eigen::MatrixXd dSubdq_plus =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_plus_ref =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus_ref =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  robot.dSubtractdConfigurationPlus(q_plus, q_minus, dSubdq_plus);
  robot.dSubtractdConfigurationMinus(q_plus, q_minus, dSubdq_minus);
  pinocchio::dDifference(model, q_minus, q_plus, dSubdq_plus_ref, pinocchio::ARG1);
  pinocchio::dDifference(model, q_minus, q_plus, dSubdq_minus_ref, pinocchio::ARG0);
  EXPECT_TRUE(dSubdq_plus.isApprox(dSubdq_plus_ref));
  EXPECT_TRUE(dSubdq_minus.isApprox(dSubdq_minus_ref));
  Eigen::MatrixXd dSubdq_plus_inv =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dSubdq_minus_inv =  Eigen::MatrixXd::Zero(model.nv, model.nv);
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationInverse(dSubdq_plus, dSubdq_plus_inv);
    robot.dSubtractdConfigurationInverse(dSubdq_minus, dSubdq_minus_inv);
    EXPECT_TRUE((dSubdq_plus.topLeftCorner(6, 6)*dSubdq_plus_inv.topLeftCorner(6, 6)).isIdentity());
    EXPECT_TRUE((dSubdq_minus.topLeftCorner(6, 6)*dSubdq_minus_inv.topLeftCorner(6, 6)).isIdentity());
    std::cout << dSubdq_plus << std::endl;
    std::cout << dSubdq_plus_inv << std::endl;
    std::cout << dSubdq_plus * dSubdq_plus_inv  << std::endl;
    std::cout << dSubdq_minus << std::endl;
    std::cout << dSubdq_minus_inv << std::endl;
    std::cout << dSubdq_minus * dSubdq_minus_inv << std::endl;
  }
}


void RobotTest::testFrameKinematics(const std::string& path_to_urdf, 
                                    const BaseJointType& base_joint_type,
                                    pinocchio::Model& model, 
                                    pinocchio::Data& data, const int frame_id) const {
  Robot robot(path_to_urdf, base_joint_type);
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


void RobotTest::testBaumgarte(const std::string& path_to_urdf, 
                              const BaseJointType& base_joint_type,
                              pinocchio::Model& model, pinocchio::Data& data, 
                              const std::vector<int>& frames) const {
  std::vector<PointContact> contacts_ref; 
  for (int i=0; i<frames.size(); ++i) {
    contacts_ref.push_back(PointContact(model, frames[i], baumgarte_weights.first, baumgarte_weights.second));
  }
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  ASSERT_EQ(robot.dimq(), model.nq);
  ASSERT_EQ(robot.dimv(), model.nv);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  robot.updateKinematics(q, v, a);
  robot.computeBaumgarteResidual(contact_status, contact_status.contactPoints(), residual.head(contact_status.dimf()));
  pinocchio::forwardKinematics(model, data, q, v, a);
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, a);
  int num_active_contacts = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      contacts_ref[i].computeBaumgarteResidual(
          model, data, contact_status.contactPoints()[i], 
          residual_ref.segment<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::MatrixXd baumgarte_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  num_active_contacts = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      contacts_ref[i].computeBaumgarteDerivatives(
          model, data, baumgarte_partial_q_ref.middleRows<3>(3*num_active_contacts), 
          baumgarte_partial_v_ref.middleRows<3>(3*num_active_contacts), 
          baumgarte_partial_a_ref.middleRows<3>(3*num_active_contacts));
      ++num_active_contacts;
    }
  }
  Eigen::MatrixXd baumgarte_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeBaumgarteDerivatives(contact_status, 
                                    baumgarte_partial_q.topRows(contact_status.dimf()), 
                                    baumgarte_partial_v.topRows(contact_status.dimf()), 
                                    baumgarte_partial_a.topRows(contact_status.dimf()));
  EXPECT_TRUE(baumgarte_partial_q.isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(baumgarte_partial_v.isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(baumgarte_partial_a.isApprox(baumgarte_partial_a_ref));
}


void RobotTest::testBaumgarteTimeStep(const std::string& path_to_urdf, 
                                      const BaseJointType& base_joint_type,
                                      pinocchio::Model& model, pinocchio::Data& data, 
                                      const std::vector<int>& frames) const {
  const double time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
  const double baumgarte_weight_velocity = 2.0 / time_step;
  const double baumgarte_weight_position = 1.0 / (time_step*time_step);
  Robot robot_time_step(path_to_urdf, base_joint_type, frames, time_step);
  Robot robot(path_to_urdf, base_joint_type, frames, 
              std::make_pair(baumgarte_weight_velocity, baumgarte_weight_position));
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(contact_status.dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(contact_status.dimf());
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  robot_time_step.updateKinematics(q, v, a);
  robot_time_step.computeBaumgarteResidual(contact_status, contact_status.contactPoints(), residual);
  robot.updateKinematics(q, v, a);
  robot.computeBaumgarteResidual(contact_status, contact_status.contactPoints(), residual_ref);
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::MatrixXd baumgarte_partial_q_ref = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v_ref = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a_ref = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_q = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_v = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  Eigen::MatrixXd baumgarte_partial_a = Eigen::MatrixXd::Zero(contact_status.dimf(), model.nv);
  robot_time_step.computeBaumgarteDerivatives(contact_status, baumgarte_partial_q, 
                                              baumgarte_partial_v, baumgarte_partial_a);
  robot.computeBaumgarteDerivatives(contact_status, baumgarte_partial_q_ref, 
                                    baumgarte_partial_v_ref, baumgarte_partial_a_ref);
  EXPECT_TRUE(baumgarte_partial_q.isApprox(baumgarte_partial_q_ref));
  EXPECT_TRUE(baumgarte_partial_v.isApprox(baumgarte_partial_v_ref));
  EXPECT_TRUE(baumgarte_partial_a.isApprox(baumgarte_partial_a_ref));
}


void RobotTest::testImpulseVelocity(const std::string& path_to_urdf, 
                                    const BaseJointType& base_joint_type,
                                    pinocchio::Model& model, pinocchio::Data& data, 
                                    const std::vector<int>& frames) const {
  std::vector<PointContact> contacts_ref; 
  for (int i=0; i<frames.size(); ++i) {
    contacts_ref.push_back(PointContact(model, frames[i], baumgarte_weights.first, baumgarte_weights.second));
  }
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  robot.updateKinematics(q, v);
  robot.computeImpulseVelocityResidual(impulse_status, residual.head(impulse_status.dimf()));
  pinocchio::forwardKinematics(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  int num_active_impulse = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      contacts_ref[i].computeContactVelocityResidual(
          model, data, residual_ref.segment<3>(3*num_active_impulse));
      ++num_active_impulse;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::MatrixXd velocity_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd velocity_partial_v_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  num_active_impulse = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      contacts_ref[i].computeContactVelocityDerivatives(
          model, data, 
          velocity_partial_q_ref.middleRows<3>(3*num_active_impulse), 
          velocity_partial_v_ref.middleRows<3>(3*num_active_impulse));
      ++num_active_impulse;
    }
  }
  Eigen::MatrixXd velocity_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  Eigen::MatrixXd velocity_partial_v = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeImpulseVelocityDerivatives(impulse_status, velocity_partial_q.topRows(impulse_status.dimf()), 
                                          velocity_partial_v.topRows(impulse_status.dimf()));
  EXPECT_TRUE(velocity_partial_q.isApprox(velocity_partial_q_ref));
  EXPECT_TRUE(velocity_partial_v.isApprox(velocity_partial_v_ref));
}


void RobotTest::testContactPosition(const std::string& path_to_urdf, 
                                    const BaseJointType& base_joint_type, 
                                    pinocchio::Model& model, pinocchio::Data& data, 
                                    const std::vector<int>& frames) const {
  std::vector<PointContact> contacts_ref; 
  for (int i=0; i<frames.size(); ++i) {
    contacts_ref.push_back(PointContact(model, frames[i], baumgarte_weights.first, baumgarte_weights.second));
  }
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(robot.max_dimf());
  Eigen::VectorXd residual_ref = Eigen::VectorXd::Zero(robot.max_dimf());
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  for (int i=0; i<frames.size(); ++i) {
    impulse_status.setContactPoint(i, Eigen::Vector3d::Random());
  }
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
  robot.updateKinematics(q, v);
  robot.computeContactPositionResidual(impulse_status, impulse_status.contactPoints(), residual.head(impulse_status.dimf()));
  pinocchio::forwardKinematics(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  pinocchio::updateFramePlacements(model, data);
  pinocchio::computeForwardKinematicsDerivatives(model, data, q, v, Eigen::VectorXd::Zero(model.nv));
  int num_active_impulse = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      contacts_ref[i].computeContactPositionResidual(
        model, data, impulse_status.contactPoints()[i], 
        residual_ref.segment<3>(3*num_active_impulse));
        ++num_active_impulse;
    }
  }
  EXPECT_TRUE(residual.isApprox(residual_ref));
  Eigen::MatrixXd contact_partial_q_ref = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  num_active_impulse = 0;
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      contacts_ref[i].computeContactPositionDerivative(
          model, data, contact_partial_q_ref.middleRows<3>(3*num_active_impulse));
      ++num_active_impulse;
    }
  }
  Eigen::MatrixXd contact_partial_q = Eigen::MatrixXd::Zero(robot.max_dimf(), model.nv);
  robot.computeContactPositionDerivative(impulse_status, contact_partial_q.topRows(impulse_status.dimf()));
  EXPECT_TRUE(contact_partial_q.isApprox(contact_partial_q_ref));
}


void RobotTest::testRNEA(const std::string& path_to_urdf, 
                         const BaseJointType& base_joint_type, 
                         pinocchio::Model& model, pinocchio::Data& data) const {
  Robot robot(path_to_urdf, base_joint_type);
  const double time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
  robot.RNEA(q, v, a, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(model, data, q, v, a);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = dRNEA_da;
  robot.RNEADerivatives(q, v, a, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::computeRNEADerivatives(model, data, q, v, a, dRNEA_dq_ref, 
                                    dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
}


void RobotTest::testRNEA(const std::string& path_to_urdf, 
                         const BaseJointType& base_joint_type, 
                         pinocchio::Model& model, pinocchio::Data& data, 
                         const std::vector<int>& frames) const {
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  std::vector<PointContact> contacts_ref; 
  for (int i=0; i<frames.size(); ++i) {
    contacts_ref.push_back(PointContact(model, frames[i], baumgarte_weights.first, baumgarte_weights.second));
  }
  const double time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd v = Eigen::VectorXd::Random(model.nv);
  const Eigen::VectorXd a = Eigen::VectorXd::Random(model.nv);
  std::vector<Eigen::Vector3d> f;
  for (const auto frame : frames) {
    f.push_back(Eigen::Vector3d::Random());
  }
  auto contact_status = robot.createContactStatus();
  contact_status.setRandom();
  robot.setContactForces(contact_status, f);
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(model.nv);
  robot.RNEA(q, v, a, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (contact_status.isContactActive(i)) {
      contacts_ref[i].computeJointForceFromContactForce(f[i], fjoint);
    }
  }
  const Eigen::VectorXd tau_ref = pinocchio::rnea(model, data, q, v, a, fjoint);
  EXPECT_TRUE(tau_ref.isApprox(tau));
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_da = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_dv;
  Eigen::MatrixXd dRNEA_da_ref = dRNEA_da;
  robot.RNEADerivatives(q, v, a, dRNEA_dq, dRNEA_dv, dRNEA_da);
  pinocchio::computeRNEADerivatives(model, data, q, v, a, fjoint, dRNEA_dq_ref, 
                                    dRNEA_dv_ref, dRNEA_da_ref);
  dRNEA_da_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_da_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_dv.isApprox(dRNEA_dv_ref));
  EXPECT_TRUE(dRNEA_da.isApprox(dRNEA_da_ref));
}


void RobotTest::testRNEAImpulse(const std::string& path_to_urdf, 
                                const BaseJointType& base_joint_type, 
                                pinocchio::Model& model, pinocchio::Data& data, 
                                const std::vector<int>& frames) const {
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  std::vector<PointContact> contacts_ref; 
  for (int i=0; i<frames.size(); ++i) {
    contacts_ref.push_back(PointContact(model, frames[i], baumgarte_weights.first, baumgarte_weights.second));
  }
  const Eigen::VectorXd q = pinocchio::randomConfiguration(
      model, -Eigen::VectorXd::Ones(model.nq), Eigen::VectorXd::Ones(model.nq));
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(model.nv);
  std::vector<Eigen::Vector3d> f;
  auto impulse_status = robot.createImpulseStatus();
  impulse_status.setRandom();
  for (const auto frame : frames) {
    f.push_back(Eigen::Vector3d::Random());
  }
  robot.setImpulseForces(impulse_status, f);
  Eigen::VectorXd impulse_res = Eigen::VectorXd::Zero(model.nv);
  robot.RNEAImpulse(q, dv, impulse_res);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model.joints.size(), pinocchio::Force::Zero());
  for (int i=0; i<contacts_ref.size(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      contacts_ref[i].computeJointForceFromContactForce(f[i], fjoint);
    }
  }
  auto model_zero_grav = model;
  model_zero_grav.gravity.linear().setZero();
  pinocchio::Data data_zero_grav(model_zero_grav);
  const Eigen::VectorXd impulse_res_ref = pinocchio::rnea(model_zero_grav, data_zero_grav, q, Eigen::VectorXd::Zero(model.nv), dv, fjoint);
  EXPECT_TRUE(impulse_res_ref.isApprox(impulse_res));
  Eigen::MatrixXd dRNEA_dq = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_ddv = Eigen::MatrixXd::Zero(model.nv, model.nv);
  Eigen::MatrixXd dRNEA_dq_ref = dRNEA_dq;
  Eigen::MatrixXd dRNEA_dv_ref = dRNEA_ddv;
  Eigen::MatrixXd dRNEA_ddv_ref = dRNEA_ddv;
  robot.RNEAImpulseDerivatives(q, dv, dRNEA_dq, dRNEA_ddv);
  pinocchio::computeRNEADerivatives(model_zero_grav, data_zero_grav, q, Eigen::VectorXd::Zero(model.nv), dv, fjoint, dRNEA_dq_ref, 
                                    dRNEA_dv_ref, dRNEA_ddv_ref);
  dRNEA_ddv_ref.triangularView<Eigen::StrictlyLower>() 
      = dRNEA_ddv_ref.transpose().triangularView<Eigen::StrictlyLower>();
  EXPECT_TRUE(dRNEA_dq.isApprox(dRNEA_dq_ref));
  EXPECT_TRUE(dRNEA_ddv.isApprox(dRNEA_ddv_ref));
}


void RobotTest::testMJtJinv(const std::string& path_to_urdf, 
                            const BaseJointType& base_joint_type, 
                            pinocchio::Model& model, pinocchio::Data& data, 
                            const std::vector<int>& frames) const {
  Robot robot(path_to_urdf, base_joint_type, frames, baumgarte_weights);
  const double time_step = std::abs(Eigen::VectorXd::Random(1)[0]);
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
  Eigen::MatrixXd dRNEA_dfext = Eigen::MatrixXd::Zero(model.nv, dimf);
  robot.updateKinematics(q, v, a);
  robot.RNEADerivatives(q, v, a, dRNEA_dq, dRNEA_dv, dRNEA_da);
  robot.dRNEAPartialdFext(contact_status, dRNEA_dfext);
  Eigen::MatrixXd MJtJinv = Eigen::MatrixXd::Zero(model.nv+dimf, model.nv+dimf);
  const Eigen::MatrixXd J = - dRNEA_dfext.transpose();
  robot.computeMJtJinv(dRNEA_da, J, MJtJinv);
  Eigen::MatrixXd MJtJ = Eigen::MatrixXd::Zero(model.nv+dimf, model.nv+dimf);
  MJtJ.topLeftCorner(model.nv, model.nv) = dRNEA_da;
  MJtJ.topRightCorner(model.nv, dimf) = J.transpose();
  MJtJ.bottomLeftCorner(dimf, model.nv) = J;
  const Eigen::MatrixXd MJtJinv_ref = MJtJ.inverse();
  EXPECT_TRUE(MJtJinv.isApprox(MJtJinv_ref));
  EXPECT_TRUE((MJtJinv*MJtJ).isIdentity());
  Eigen::MatrixXd Minv = dRNEA_da;
  robot.computeMinv(dRNEA_da, Minv);
  EXPECT_TRUE((dRNEA_da*Minv).isIdentity());
  EXPECT_TRUE((Minv*dRNEA_da).isIdentity());
  std::cout << Minv << std::endl;
}


void RobotTest::testGenConfiguration(const std::string& path_to_urdf, 
                                     const BaseJointType& base_joint_type, 
                                     pinocchio::Model& model, 
                                     pinocchio::Data& data) const {
  Robot robot(path_to_urdf, base_joint_type);
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


TEST_F(RobotTest, hasFloatingBase) {
  Robot fixed_base(fixed_base_urdf);
  EXPECT_FALSE(fixed_base.hasFloatingBase());
  EXPECT_EQ(fixed_base.dim_passive(), 0);
  EXPECT_EQ(fixed_base.dimu(), fixed_base.dimv());
  Robot fixed_base_contacts(fixed_base_urdf, BaseJointType::FixedBase, 
                            fixed_base_contact_frames, baumgarte_weights);
  EXPECT_FALSE(fixed_base_contacts.hasFloatingBase());
  EXPECT_EQ(fixed_base_contacts.dim_passive(), 0);
  EXPECT_EQ(fixed_base_contacts.dimu(), fixed_base_contacts.dimv());
  Robot floating_base(floating_base_urdf, BaseJointType::FloatingBase);
  EXPECT_TRUE(floating_base.hasFloatingBase());
  EXPECT_EQ(floating_base.dim_passive(), 6);
  EXPECT_EQ(floating_base.dimu()+floating_base.dim_passive(), floating_base.dimv());
  Robot floating_base_contacts(floating_base_urdf, BaseJointType::FloatingBase, 
                               floating_base_contact_frames, baumgarte_weights);
  EXPECT_TRUE(floating_base_contacts.hasFloatingBase());
  EXPECT_EQ(floating_base_contacts.dim_passive(), 6);
  EXPECT_EQ(floating_base_contacts.dimu()+floating_base_contacts.dim_passive(), 
            floating_base_contacts.dimv());
}


TEST_F(RobotTest, testFixedbase) {
  const auto path_to_urdf = fixed_base_urdf;
  const auto contact_frames = fixed_base_contact_frames;
  pinocchio::Model model = fixed_base_robot;
  pinocchio::Data data = fixed_base_data;
  testConstructorAndSetter(path_to_urdf, BaseJointType::FixedBase, contact_frames);
  testIntegrateConfiguration(path_to_urdf, BaseJointType::FixedBase, model);
  testSubtractConfiguration(path_to_urdf, BaseJointType::FixedBase, model);
  testSubtractConfigurationDerivatives(path_to_urdf, BaseJointType::FixedBase, model);
  for (const auto frame : contact_frames) {
    testFrameKinematics(path_to_urdf, BaseJointType::FixedBase, model, data, frame);
  }
  testBaumgarte(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testImpulseVelocity(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testContactPosition(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testRNEA(path_to_urdf, BaseJointType::FixedBase, model, data);
  testRNEA(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testRNEAImpulse(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testMJtJinv(path_to_urdf, BaseJointType::FixedBase, model, data, contact_frames);
  testGenConfiguration(path_to_urdf, BaseJointType::FixedBase, model, data);
}


TEST_F(RobotTest, testFloatingBase) {
  const auto path_to_urdf = floating_base_urdf;
  const auto contact_frames = floating_base_contact_frames;
  pinocchio::Model model = floating_base_robot;
  pinocchio::Data data = floating_base_data;
  testConstructorAndSetter(path_to_urdf, BaseJointType::FloatingBase, contact_frames);
  testIntegrateConfiguration(path_to_urdf, BaseJointType::FloatingBase, model);
  testSubtractConfiguration(path_to_urdf, BaseJointType::FloatingBase, model);
  testSubtractConfigurationDerivatives(path_to_urdf, BaseJointType::FloatingBase, model);
  for (const auto frame : contact_frames) {
    testFrameKinematics(path_to_urdf, BaseJointType::FloatingBase, model, data, frame);
  }
  testBaumgarte(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testImpulseVelocity(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testContactPosition(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testRNEA(path_to_urdf, BaseJointType::FloatingBase, model, data);
  testRNEA(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testRNEAImpulse(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testMJtJinv(path_to_urdf, BaseJointType::FloatingBase, model, data, contact_frames);
  testGenConfiguration(path_to_urdf, BaseJointType::FloatingBase, model, data);
}

} // namespace idocp 


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}