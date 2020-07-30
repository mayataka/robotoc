#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "constraints/joint_space_constraints/joint_position_lower_limits.hpp"
#include "constraints/pdipm_func.hpp"
#include "robot/robot.hpp"


namespace idocp {
namespace pdipm {


class JointPositionLowerLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    const std::string urdf_file_name = "../../../urdf/iiwa14/iiwa14.urdf";
    robot_ = Robot(urdf_file_name);
    dimq_ = robot_.dimq();
    dimv_ = robot_.dimv();
    barrier_ = 1.0e-04;
    dtau_ = 0.1;
    qmin_ = robot_.lowerJointPositionLimit();
    slack_.array() = Eigen::VectorXd::Random(dimq_).array().abs();
    dual_.array() = Eigen::VectorXd::Random(dimq_).array().abs();
    dslack_.array() = Eigen::VectorXd::Random(dimq_);
    ddual_.array() = Eigen::VectorXd::Random(dimq_);
  }

  virtual void TearDown() {
  }

  Robot robot_;
  int dimq_, dimv_;
  double barrier_, dtau_;
  Eigen::VectorXd qmin_, slack_, dual_, dslack_, ddual_;
};


TEST_F(JointPositionLowerLimitTest, isFeasible) {
  JointPositionLowerLimits limit_(robot_, barrier_);
  Eigen::VectorXd q = qmin_;
  const Eigen::VectorXd tmp = Eigen::VectorXd::Random(dimq_);
  q += tmp;
  if (tmp.minCoeff() < 0) {
    EXPECT_FALSE(limit_.isFeasible(robot_, q));
  }
  else {
    EXPECT_TRUE(limit_.isFeasible(robot_, q));
  }
}


TEST_F(JointPositionLowerLimitTest, setSlackAndDual) {
  JointPositionLowerLimits limit_(robot_, barrier_);
  Eigen::VectorXd tmp = Eigen::VectorXd::Random(dimq_);
  const Eigen::VectorXd q = qmin_ + tmp;
  limit_.setSlackAndDual(robot_, dtau_, q);
  Eigen::VectorXd slack_ref = dtau_ * (q-qmin_);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq_);
  pdipmfunc::SetSlackAndDualPositive(dimq_, barrier_, slack_ref, dual_ref);
  EXPECT_DOUBLE_EQ(limit_.slackBarrier(), 
                   pdipmfunc::SlackBarrierCost(dimq_, barrier_, slack_ref));
  Eigen::VectorXd dual = Eigen::VectorXd::Zero(dimq_);
  limit_.augmentDualResidual(robot_, dtau_, dual);
  dual_ref *= dtau_;
  EXPECT_TRUE(dual.isApprox(dual_ref));
}


TEST_F(JointPositionLowerLimitTest, condenseSlackAndDual) {
  JointPositionLowerLimits limit_(robot_, barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Random(dimq_);
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::VectorXd Cq = Eigen::VectorXd::Zero(dimq_);
  limit_.setSlackAndDual(robot_, dtau_, q);
  limit_.condenseSlackAndDual(robot_, dtau_, q, Cqq, Cq);
  Eigen::VectorXd slack_ref = dtau_ * (q-qmin_);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq_);
  pdipmfunc::SetSlackAndDualPositive(dimq_, barrier_, slack_ref, dual_ref);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Zero(dimq_);
  for (int i=0; i<dimv_; ++i) {
    Cqq_ref(i, i) = dtau_ * dtau_ * dual_ref(i) / slack_ref(i);
  }
  Eigen::VectorXd residual = dtau_ * (qmin_-q) + slack_ref;
  Eigen::VectorXd duality = Eigen::VectorXd::Zero(dimq_);
  for (int i=0; i<dimq_; ++i) {
    duality(i) = dual_ref(i) * slack_ref(i) - barrier_;
  }
  for (int i=0; i<dimq_; ++i) {
    Cq_ref(i) = - dtau_ * (dual_ref(i)*residual(i)-duality(i)) / slack_ref(i);
  }
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
}


TEST_F(JointPositionLowerLimitTest, computeDirectionAndMaxStepSize) {
  JointPositionLowerLimits limit_(robot_, barrier_);
  Eigen::VectorXd q = Eigen::VectorXd::Random(dimq_);
  Eigen::MatrixXd Cqq = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::VectorXd Cq = Eigen::VectorXd::Zero(dimq_);
  limit_.setSlackAndDual(robot_, dtau_, q);
  limit_.condenseSlackAndDual(robot_, dtau_, q, Cqq, Cq);
  Eigen::VectorXd slack_ref = dtau_ * (q-qmin_);
  Eigen::VectorXd dual_ref = Eigen::VectorXd::Zero(dimq_);
  pdipmfunc::SetSlackAndDualPositive(dimq_, barrier_, slack_ref, dual_ref);
  Eigen::MatrixXd Cqq_ref = Eigen::MatrixXd::Zero(dimq_, dimq_);
  Eigen::VectorXd Cq_ref = Eigen::VectorXd::Zero(dimq_);
  for (int i=0; i<dimv_; ++i) {
    Cqq_ref(i, i) = dtau_ * dtau_ * dual_ref(i) / slack_ref(i);
  }
  Eigen::VectorXd residual = dtau_ * (qmin_-q) + slack_ref;
  Eigen::VectorXd duality = Eigen::VectorXd::Zero(dimq_);
  for (int i=0; i<dimq_; ++i) {
    duality(i) = dual_ref(i) * slack_ref(i) - barrier_;
  }
  for (int i=0; i<dimq_; ++i) {
    Cq_ref(i) = - dtau_ * (dual_ref(i)*residual(i)-duality(i)) / slack_ref(i);
  }
  EXPECT_TRUE(Cqq.isApprox(Cqq_ref));
  EXPECT_TRUE(Cq.isApprox(Cq_ref));
}



} // namespace pdipm
} // namespace invdynocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}