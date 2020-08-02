#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/joint_space_constraints/joint_space_constraints.hpp"
#include "constraints/joint_space_constraints/joint_variables_lower_limits.hpp"
#include "constraints/joint_space_constraints/joint_variables_upper_limits.hpp"


namespace idocp {
namespace pdipm {

class JointSpaceConstraintsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf_ = "../../../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf_ = "../../../urdf/anymal/anymal.urdf";
    fixed_base_robot_ = Robot(fixed_base_urdf_);
    floating_base_robot_ = Robot(floating_base_urdf_);
    barrier_ = 1.0e-04;
    dtau_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    margin_rate_ = 0.995;
    while (barrier_ == 0) {
      barrier_ = std::abs(Eigen::VectorXd::Random(1)[0]);
    }
  }

  virtual void TearDown() {
  }

  double barrier_, dtau_, margin_rate_;
  Eigen::VectorXd slack_, dual_, dslack_, ddual_;
  std::string fixed_base_urdf_, floating_base_urdf_;
  Robot fixed_base_robot_, floating_base_robot_;
};


TEST_F(JointSpaceConstraintsTest, isFeasibleFixedBase) {
  Robot robot = fixed_base_robot_;
  JointSpaceConstraints constraints(robot);
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  JointVariablesUpperLimits position_upper_limits_ref(robot, robot.upperJointPositionLimit(), barrier_);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits position_lower_limits_ref(robot, robot.lowerJointPositionLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  if (!position_upper_limits_ref.isFeasible(q)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!velocity_upper_limits_ref.isFeasible(v)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!torque_upper_limits_ref.isFeasible(u)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!position_lower_limits_ref.isFeasible(q)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!velocity_lower_limits_ref.isFeasible(v)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!torque_lower_limits_ref.isFeasible(u)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else {
    EXPECT_TRUE(constraints.isFeasible(q, v, a, u));
  }
}


TEST_F(JointSpaceConstraintsTest, isFeasibleFloatingBase) {
  Robot robot = floating_base_robot_;
  JointSpaceConstraints constraints(robot);
  Eigen::VectorXd q = Eigen::VectorXd::Random(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  JointVariablesUpperLimits position_upper_limits_ref(robot, robot.upperJointPositionLimit(), barrier_);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits position_lower_limits_ref(robot, robot.lowerJointPositionLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  if (!position_upper_limits_ref.isFeasible(q)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!velocity_upper_limits_ref.isFeasible(v)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!torque_upper_limits_ref.isFeasible(u)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!position_lower_limits_ref.isFeasible(q)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!velocity_lower_limits_ref.isFeasible(v)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else if (!torque_lower_limits_ref.isFeasible(u)) {
    EXPECT_FALSE(constraints.isFeasible(q, v, a, u));
  }
  else {
    EXPECT_TRUE(constraints.isFeasible(q, v, a, u));
  }
}


TEST_F(JointSpaceConstraintsTest, condensingFixedBase) {
  Robot robot = fixed_base_robot_;
  JointSpaceConstraints constraints(robot);
  JointVariablesUpperLimits position_upper_limits_ref(robot, robot.upperJointPositionLimit(), barrier_);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits position_lower_limits_ref(robot, robot.lowerJointPositionLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  while (!velocity_upper_limits_ref.isFeasible(v) || !velocity_lower_limits_ref.isFeasible(v)) {
    v = Eigen::VectorXd::Random(robot.dimv());
  }
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  while (!torque_upper_limits_ref.isFeasible(u) || !torque_lower_limits_ref.isFeasible(u)) {
    u = Eigen::VectorXd::Random(robot.dimv());
  }
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  position_upper_limits_ref.setSlackAndDual(dtau_, q);
  position_lower_limits_ref.setSlackAndDual(dtau_, q);
  velocity_upper_limits_ref.setSlackAndDual(dtau_, v);
  velocity_lower_limits_ref.setSlackAndDual(dtau_, v);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  position_upper_limits_ref.augmentDualResidual(dtau_, Cq_ref);
  position_lower_limits_ref.augmentDualResidual(dtau_, Cq_ref);
  velocity_upper_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  velocity_lower_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  position_upper_limits_ref.condenseSlackAndDual(dtau_, q, Cqq_ref, Cq_ref);
  position_lower_limits_ref.condenseSlackAndDual(dtau_, q, Cqq_ref, Cq_ref);
  velocity_upper_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  velocity_lower_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  position_upper_limits_ref.computeSlackAndDualDirection(dtau_, dq);
  position_lower_limits_ref.computeSlackAndDualDirection(dtau_, dq);
  velocity_upper_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  velocity_lower_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_position_upper_limit 
      = position_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_position_lower_limit 
      = position_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_position_upper_limit, 
                                              slack_size_position_lower_limit, 
                                              slack_size_velocity_upper_limit, 
                                              slack_size_velocity_lower_limit,
                                              slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_position_upper_limit 
      = position_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_position_lower_limit 
      = position_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_position_upper_limit, 
                                             dual_size_position_lower_limit, 
                                             dual_size_velocity_upper_limit, 
                                             dual_size_velocity_lower_limit,
                                             dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier()
                                  + position_lower_limits_ref.costSlackBarrier()
                                  + velocity_upper_limits_ref.costSlackBarrier()
                                  + velocity_lower_limits_ref.costSlackBarrier()
                                  + torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + position_lower_limits_ref.costSlackBarrier(max_step_size)
                                  + velocity_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + velocity_lower_limits_ref.costSlackBarrier(max_step_size)
                                  + torque_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  position_upper_limits_ref.updateSlack(max_step_size);
  position_lower_limits_ref.updateSlack(max_step_size);
  velocity_upper_limits_ref.updateSlack(max_step_size);
  velocity_lower_limits_ref.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier()
                            + position_lower_limits_ref.costSlackBarrier()
                            + velocity_upper_limits_ref.costSlackBarrier()
                            + velocity_lower_limits_ref.costSlackBarrier()
                            + torque_upper_limits_ref.costSlackBarrier()
                            + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  position_upper_limits_ref.updateDual(max_step_size);
  position_lower_limits_ref.updateDual(max_step_size);
  velocity_upper_limits_ref.updateDual(max_step_size);
  velocity_lower_limits_ref.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = position_upper_limits_ref.residualL1Nrom(dtau_, q)
                            + position_lower_limits_ref.residualL1Nrom(dtau_, q)
                            + velocity_upper_limits_ref.residualL1Nrom(dtau_, v)
                            + velocity_lower_limits_ref.residualL1Nrom(dtau_, v)
                            + torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = position_upper_limits_ref.residualSquaredNrom(dtau_, q)
                            + position_lower_limits_ref.residualSquaredNrom(dtau_, q)
                            + velocity_upper_limits_ref.residualSquaredNrom(dtau_, v)
                            + velocity_lower_limits_ref.residualSquaredNrom(dtau_, v)
                            + torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}


TEST_F(JointSpaceConstraintsTest, condensingFloatingBase) {
  Robot robot = floating_base_robot_;
  JointSpaceConstraints constraints(robot);
  JointVariablesUpperLimits position_upper_limits_ref(robot, robot.upperJointPositionLimit(), barrier_);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits position_lower_limits_ref(robot, robot.lowerJointPositionLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  while (!velocity_upper_limits_ref.isFeasible(v) || !velocity_lower_limits_ref.isFeasible(v)) {
    v = Eigen::VectorXd::Random(robot.dimv());
  }
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  while (!torque_upper_limits_ref.isFeasible(u) || !torque_lower_limits_ref.isFeasible(u)) {
    u = Eigen::VectorXd::Random(robot.dimv());
  }
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  position_upper_limits_ref.setSlackAndDual(dtau_, q);
  position_lower_limits_ref.setSlackAndDual(dtau_, q);
  velocity_upper_limits_ref.setSlackAndDual(dtau_, v);
  velocity_lower_limits_ref.setSlackAndDual(dtau_, v);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  position_upper_limits_ref.augmentDualResidual(dtau_, Cq_ref);
  position_lower_limits_ref.augmentDualResidual(dtau_, Cq_ref);
  velocity_upper_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  velocity_lower_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  position_upper_limits_ref.condenseSlackAndDual(dtau_, q, Cqq_ref, Cq_ref);
  position_lower_limits_ref.condenseSlackAndDual(dtau_, q, Cqq_ref, Cq_ref);
  velocity_upper_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  velocity_lower_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  position_upper_limits_ref.computeSlackAndDualDirection(dtau_, dq);
  position_lower_limits_ref.computeSlackAndDualDirection(dtau_, dq);
  velocity_upper_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  velocity_lower_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_position_upper_limit 
      = position_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_position_lower_limit 
      = position_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_position_upper_limit, 
                                              slack_size_position_lower_limit, 
                                              slack_size_velocity_upper_limit, 
                                              slack_size_velocity_lower_limit,
                                              slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_position_upper_limit 
      = position_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_position_lower_limit 
      = position_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_position_upper_limit, 
                                             dual_size_position_lower_limit, 
                                             dual_size_velocity_upper_limit, 
                                             dual_size_velocity_lower_limit,
                                             dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier()
                                  + position_lower_limits_ref.costSlackBarrier()
                                  + velocity_upper_limits_ref.costSlackBarrier()
                                  + velocity_lower_limits_ref.costSlackBarrier()
                                  + torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + position_lower_limits_ref.costSlackBarrier(max_step_size)
                                  + velocity_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + velocity_lower_limits_ref.costSlackBarrier(max_step_size)
                                  + torque_upper_limits_ref.costSlackBarrier(max_step_size)
                                  + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  position_upper_limits_ref.updateSlack(max_step_size);
  position_lower_limits_ref.updateSlack(max_step_size);
  velocity_upper_limits_ref.updateSlack(max_step_size);
  velocity_lower_limits_ref.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = position_upper_limits_ref.costSlackBarrier()
                            + position_lower_limits_ref.costSlackBarrier()
                            + velocity_upper_limits_ref.costSlackBarrier()
                            + velocity_lower_limits_ref.costSlackBarrier()
                            + torque_upper_limits_ref.costSlackBarrier()
                            + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  position_upper_limits_ref.updateDual(max_step_size);
  position_lower_limits_ref.updateDual(max_step_size);
  velocity_upper_limits_ref.updateDual(max_step_size);
  velocity_lower_limits_ref.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = position_upper_limits_ref.residualL1Nrom(dtau_, q)
                            + position_lower_limits_ref.residualL1Nrom(dtau_, q)
                            + velocity_upper_limits_ref.residualL1Nrom(dtau_, v)
                            + velocity_lower_limits_ref.residualL1Nrom(dtau_, v)
                            + torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = position_upper_limits_ref.residualSquaredNrom(dtau_, q)
                            + position_lower_limits_ref.residualSquaredNrom(dtau_, q)
                            + velocity_upper_limits_ref.residualSquaredNrom(dtau_, v)
                            + velocity_lower_limits_ref.residualSquaredNrom(dtau_, v)
                            + torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}


TEST_F(JointSpaceConstraintsTest, condensingFixedBaseTimeStep0) {
  Robot robot = fixed_base_robot_;
  JointSpaceConstraints constraints(robot);
  constraints.setTimeStep(0);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier(max_step_size)
                           + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier()
                           + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}


TEST_F(JointSpaceConstraintsTest, condensingFloatingBaseTimeStep0) {
  Robot robot = floating_base_robot_;
  JointSpaceConstraints constraints(robot);
  constraints.setTimeStep(0);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier(max_step_size)
                           + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = torque_upper_limits_ref.costSlackBarrier()
                           + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}


TEST_F(JointSpaceConstraintsTest, condensingFixedBaseTimeStep1) {
  Robot robot = fixed_base_robot_;
  JointSpaceConstraints constraints(robot);
  constraints.setTimeStep(1);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  while (!velocity_upper_limits_ref.isFeasible(v) || !velocity_lower_limits_ref.isFeasible(v)) {
    v = Eigen::VectorXd::Random(robot.dimv());
  }
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  while (!torque_upper_limits_ref.isFeasible(u) || !torque_lower_limits_ref.isFeasible(u)) {
    u = Eigen::VectorXd::Random(robot.dimv());
  }
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  velocity_upper_limits_ref.setSlackAndDual(dtau_, v);
  velocity_lower_limits_ref.setSlackAndDual(dtau_, v);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  velocity_upper_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  velocity_lower_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  velocity_upper_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  velocity_lower_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  velocity_upper_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  velocity_lower_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_velocity_upper_limit, 
                                              slack_size_velocity_lower_limit,
                                              slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_velocity_upper_limit, 
                                             dual_size_velocity_lower_limit,
                                             dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier()
                                  + velocity_lower_limits_ref.costSlackBarrier()
                                  + torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier(max_step_size)
                           + velocity_lower_limits_ref.costSlackBarrier(max_step_size)
                           + torque_upper_limits_ref.costSlackBarrier(max_step_size)
                           + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  velocity_upper_limits_ref.updateSlack(max_step_size);
  velocity_lower_limits_ref.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier()
                           + velocity_lower_limits_ref.costSlackBarrier()
                           + torque_upper_limits_ref.costSlackBarrier()
                           + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  velocity_upper_limits_ref.updateDual(max_step_size);
  velocity_lower_limits_ref.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = velocity_upper_limits_ref.residualL1Nrom(dtau_, v)
                            + velocity_lower_limits_ref.residualL1Nrom(dtau_, v)
                            + torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = velocity_upper_limits_ref.residualSquaredNrom(dtau_, v)
                            + velocity_lower_limits_ref.residualSquaredNrom(dtau_, v)
                            + torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}


TEST_F(JointSpaceConstraintsTest, condensingFloatingBaseTimeStep1) {
  Robot robot = floating_base_robot_;
  JointSpaceConstraints constraints(robot);
  constraints.setTimeStep(1);
  JointVariablesUpperLimits velocity_upper_limits_ref(robot, robot.jointVelocityLimit(), barrier_);
  JointVariablesUpperLimits torque_upper_limits_ref(robot, robot.jointEffortLimit(), barrier_);
  JointVariablesLowerLimits velocity_lower_limits_ref(robot, -robot.jointVelocityLimit(), barrier_);
  JointVariablesLowerLimits torque_lower_limits_ref(robot, -robot.jointEffortLimit(), barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(robot.dimq());
  robot.generateFeasibleConfiguration(q);
  Eigen::VectorXd v = Eigen::VectorXd::Random(robot.dimv());
  while (!velocity_upper_limits_ref.isFeasible(v) || !velocity_lower_limits_ref.isFeasible(v)) {
    v = Eigen::VectorXd::Random(robot.dimv());
  }
  Eigen::VectorXd a = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd u = Eigen::VectorXd::Random(robot.dimv());
  while (!torque_upper_limits_ref.isFeasible(u) || !torque_lower_limits_ref.isFeasible(u)) {
    u = Eigen::VectorXd::Random(robot.dimv());
  }
  ASSERT_TRUE(constraints.isFeasible(q, v, a, u));
  constraints.setSlackAndDual(dtau_, q, v, a, u);
  velocity_upper_limits_ref.setSlackAndDual(dtau_, v);
  velocity_lower_limits_ref.setSlackAndDual(dtau_, v);
  torque_upper_limits_ref.setSlackAndDual(dtau_, u);
  torque_lower_limits_ref.setSlackAndDual(dtau_, u);
  Eigen::VectorXd Cq = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  constraints.augmentDualResidual(dtau_, Cu);
  constraints.augmentDualResidual(dtau_, Cq, Cv, Ca);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cv_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Cu_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  Eigen::VectorXd Ca_ref = 1e-05 * Eigen::VectorXd::Ones(robot.dimv());
  velocity_upper_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  velocity_lower_limits_ref.augmentDualResidual(dtau_, Cv_ref);
  torque_upper_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  torque_lower_limits_ref.augmentDualResidual(dtau_, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  std::cout << "Cq " << std::endl;
  std::cout << Cq.transpose() << "\n" << std::endl;
  std::cout << "Cq_ref " << std::endl;
  std::cout << Cq_ref.transpose() << "\n" << std::endl;
  std::cout << "Cv " << std::endl;
  std::cout << Cv.transpose() << "\n" << std::endl;
  std::cout << "Cv_ref " << std::endl;
  std::cout << Cv_ref.transpose() << "\n" << std::endl;
  std::cout << "Ca " << std::endl;
  std::cout << Ca.transpose() << "\n" << std::endl;
  std::cout << "Ca_ref " << std::endl;
  std::cout << Ca_ref.transpose() << "\n" << std::endl;
  std::cout << "Cu " << std::endl;
  std::cout << Cu.transpose() << "\n" << std::endl;
  std::cout << "Cu_ref " << std::endl;
  std::cout << Cu_ref.transpose() << "\n" << std::endl;
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  constraints.condenseSlackAndDual(dtau_, q, v, a, Cqq, Cvv, Caa, Cq, Cv, Ca);
  constraints.condenseSlackAndDual(dtau_, u, Cuu, Cu);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cvv_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Caa_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  Eigen::MatrixXd Cuu_ref = 1e-05 * Eigen::MatrixXd::Ones(robot.dimv(), robot.dimv());
  velocity_upper_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  velocity_lower_limits_ref.condenseSlackAndDual(dtau_, v, Cvv_ref, Cv_ref);
  torque_upper_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  torque_lower_limits_ref.condenseSlackAndDual(dtau_, u, Cuu_ref, Cu_ref);
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
  EXPECT_TRUE(Cv.isApprox(Cv_ref));
  EXPECT_TRUE(Ca.isApprox(Ca_ref));
  EXPECT_TRUE(Cu.isApprox(Cu_ref));
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cvv.isApprox(Cvv_ref));
  EXPECT_TRUE(Caa.isApprox(Caa_ref));
  EXPECT_TRUE(Cuu.isApprox(Cuu_ref));
  std::cout << "Cqq " << std::endl;
  std::cout << Cqq << "\n" << std::endl;
  std::cout << "Cqq_ref " << std::endl;
  std::cout << Cqq_ref << "\n" << std::endl;
  std::cout << "Cvv " << std::endl;
  std::cout << Cvv << "\n" << std::endl;
  std::cout << "Cvv_ref " << std::endl;
  std::cout << Cvv_ref << "\n" << std::endl;
  std::cout << "Caa " << std::endl;
  std::cout << Caa << "\n" << std::endl;
  std::cout << "Caa_ref " << std::endl;
  std::cout << Caa_ref << "\n" << std::endl;
  std::cout << "Cuu " << std::endl;
  std::cout << Cuu << "\n" << std::endl;
  std::cout << "Cuu_ref " << std::endl;
  std::cout << Cuu_ref << "\n" << std::endl;
  const Eigen::VectorXd dq = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd dv = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd da = Eigen::VectorXd::Random(robot.dimv());
  const Eigen::VectorXd du = Eigen::VectorXd::Random(robot.dimv());
  constraints.computeSlackAndDualDirection(dtau_, dq, dv, da, du);
  velocity_upper_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  velocity_lower_limits_ref.computeSlackAndDualDirection(dtau_, dv);
  torque_upper_limits_ref.computeSlackAndDualDirection(dtau_, du);
  torque_lower_limits_ref.computeSlackAndDualDirection(dtau_, du);
  const double max_slack_step = constraints.maxSlackStepSize();
  const double slack_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_upper_limit 
      = torque_upper_limits_ref.maxSlackStepSize(margin_rate_);
  const double slack_size_torque_lower_limit 
      = torque_lower_limits_ref.maxSlackStepSize(margin_rate_);
  const double max_slack_step_ref = std::min({slack_size_velocity_upper_limit, 
                                              slack_size_velocity_lower_limit,
                                              slack_size_torque_upper_limit, 
                                              slack_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_slack_step, max_slack_step_ref);
  const double max_dual_step = constraints.maxDualStepSize();
  const double dual_size_velocity_upper_limit 
      = velocity_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_velocity_lower_limit 
      = velocity_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_upper_limit 
      = torque_upper_limits_ref.maxDualStepSize(margin_rate_);
  const double dual_size_torque_lower_limit 
      = torque_lower_limits_ref.maxDualStepSize(margin_rate_);
  const double max_dual_step_ref = std::min({dual_size_velocity_upper_limit, 
                                             dual_size_velocity_lower_limit,
                                             dual_size_torque_upper_limit, 
                                             dual_size_torque_lower_limit});
  EXPECT_DOUBLE_EQ(max_dual_step, max_dual_step_ref);
  double cost_slack_barrier = constraints.costSlackBarrier();
  double cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier()
                                  + velocity_lower_limits_ref.costSlackBarrier()
                                  + torque_upper_limits_ref.costSlackBarrier()
                                  + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  const double max_step_size = std::min(max_slack_step, max_dual_step);
  cost_slack_barrier = constraints.costSlackBarrier(max_step_size);
  cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier(max_step_size)
                           + velocity_lower_limits_ref.costSlackBarrier(max_step_size)
                           + torque_upper_limits_ref.costSlackBarrier(max_step_size)
                           + torque_lower_limits_ref.costSlackBarrier(max_step_size);
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateSlack(max_step_size);
  velocity_upper_limits_ref.updateSlack(max_step_size);
  velocity_lower_limits_ref.updateSlack(max_step_size);
  torque_upper_limits_ref.updateSlack(max_step_size);
  torque_lower_limits_ref.updateSlack(max_step_size);
  cost_slack_barrier = constraints.costSlackBarrier();
  cost_slack_barrier_ref = velocity_upper_limits_ref.costSlackBarrier()
                           + velocity_lower_limits_ref.costSlackBarrier()
                           + torque_upper_limits_ref.costSlackBarrier()
                           + torque_lower_limits_ref.costSlackBarrier();
  EXPECT_DOUBLE_EQ(cost_slack_barrier, cost_slack_barrier_ref);
  constraints.updateDual(max_step_size);
  velocity_upper_limits_ref.updateDual(max_step_size);
  velocity_lower_limits_ref.updateDual(max_step_size);
  torque_upper_limits_ref.updateDual(max_step_size);
  torque_lower_limits_ref.updateDual(max_step_size);
  const double l1norm = constraints.residualL1Nrom(dtau_, q, v, a, u);
  const double l1norm_ref = velocity_upper_limits_ref.residualL1Nrom(dtau_, v)
                            + velocity_lower_limits_ref.residualL1Nrom(dtau_, v)
                            + torque_upper_limits_ref.residualL1Nrom(dtau_, u)
                            + torque_lower_limits_ref.residualL1Nrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l1norm, l1norm_ref);
  const double l2norm = constraints.residualSquaredNrom(dtau_, q, v, a, u);
  const double l2norm_ref = velocity_upper_limits_ref.residualSquaredNrom(dtau_, v)
                            + velocity_lower_limits_ref.residualSquaredNrom(dtau_, v)
                            + torque_upper_limits_ref.residualSquaredNrom(dtau_, u)
                            + torque_lower_limits_ref.residualSquaredNrom(dtau_, u);
  EXPECT_DOUBLE_EQ(l2norm, l2norm_ref);
}



} // namespace pdipm
} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}