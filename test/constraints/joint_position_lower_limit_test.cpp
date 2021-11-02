#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class JointPositionLowerLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void test_kinematics(Robot& robot) const;
  void test_isFeasible(Robot& robot) const;
  void test_setSlack(Robot& robot) const;
  void test_evalConstraint(Robot& robot) const;
  void test_evalDerivatives(Robot& robot) const;
  void test_condenseSlackAndDual(Robot& robot) const;
  void test_expandSlackAndDual(Robot& robot) const;

  double barrier, dt;
};


void JointPositionLowerLimitTest::test_kinematics(Robot& robot) const {
  JointPositionLowerLimit constr(robot); 
  EXPECT_FALSE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::PositionLevel);
}


void JointPositionLowerLimitTest::test_isFeasible(Robot& robot) const {
  JointPositionLowerLimit constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, data, s));
  s.q = 2*robot.lowerJointPositionLimit();
  EXPECT_FALSE(constr.isFeasible(robot, data, s));
}


void JointPositionLowerLimitTest::test_setSlack(Robot& robot) const {
  JointPositionLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter()), data_ref(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  constr.setSlack(robot, data, s);
  data_ref.slack = -qmin + s.q.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointPositionLowerLimitTest::test_evalConstraint(Robot& robot) const {
  JointPositionLowerLimit constr(robot); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  auto data_ref = data;
  constr.evalConstraint(robot, data, s);
  data_ref.residual = - s.q.tail(dimc) + qmin + data_ref.slack;
  pdipm::computeComplementarySlackness(barrier, data_ref);
  data_ref.log_barrier = pdipm::logBarrier(barrier, data_ref.slack);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointPositionLowerLimitTest::test_evalDerivatives(Robot& robot) const {
  JointPositionLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, data, s, kkt_res);
  kkt_res_ref.lq().tail(dimc) -= data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointPositionLowerLimitTest::test_condenseSlackAndDual(Robot& robot) const {
  JointPositionLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(data, s, kkt_mat, kkt_res);
  kkt_res_ref.lq().tail(dimc).array() 
      -= (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qqq().diagonal().tail(dimc).array() 
      += data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointPositionLowerLimitTest::test_expandSlackAndDual(Robot& robot) const {
  JointPositionLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmin = robot.lowerJointPositionLimit();
  constr.setSlack(robot, data, s);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  constr.expandSlackAndDual(data, s, d);
  data_ref.dslack = d.dq().tail(dimc) - data_ref.residual;
  pdipm::computeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointPositionLowerLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  test_kinematics(robot);
  test_isFeasible(robot);
  test_setSlack(robot);
  test_evalConstraint(robot);
  test_evalDerivatives(robot);
  test_condenseSlackAndDual(robot);
  test_expandSlackAndDual(robot);
}


TEST_F(JointPositionLowerLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  test_kinematics(robot);
  test_isFeasible(robot);
  test_setSlack(robot);
  test_evalConstraint(robot);
  test_evalDerivatives(robot);
  test_condenseSlackAndDual(robot);
  test_expandSlackAndDual(robot);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}