#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_position_lower_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace idocp {

class JointPositionLowerLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot) const;
  void testIsFeasible(Robot& robot) const;
  void testSetSlack(Robot& robot) const;
  void test_evalConstraint(Robot& robot) const;
  void test_evalDerivatives(Robot& robot) const;
  void testCondenseSlackAndDual(Robot& robot) const;
  void testExpandSlackAndDual(Robot& robot) const;

  double barrier, dt;
};


void JointPositionLowerLimitTest::testKinematics(Robot& robot) const {
  JointPositionLowerLimit constr(robot); 
  EXPECT_FALSE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::PositionLevel);
}


void JointPositionLowerLimitTest::testIsFeasible(Robot& robot) const {
  JointPositionLowerLimit constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, data, s));
  s.q = 2*robot.lowerJointPositionLimit();
  EXPECT_FALSE(constr.isFeasible(robot, data, s));
}


void JointPositionLowerLimitTest::testSetSlack(Robot& robot) const {
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
  pdipm::ComputeComplementarySlackness(barrier, data_ref);
  data_ref.log_barrier = pdipm::LogBarrier(barrier, data_ref.slack);
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
  constr.evalDerivatives(robot, data, dt, s, kkt_res);
  kkt_res_ref.lq().tail(dimc) -= dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointPositionLowerLimitTest::testCondenseSlackAndDual(Robot& robot) const {
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
  constr.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  kkt_res_ref.lq().tail(dimc).array() 
      -= dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qqq().diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointPositionLowerLimitTest::testExpandSlackAndDual(Robot& robot) const {
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
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointPositionLowerLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlack(robot);
  test_evalConstraint(robot);
  test_evalDerivatives(robot);
  testCondenseSlackAndDual(robot);
  testExpandSlackAndDual(robot);
}


TEST_F(JointPositionLowerLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlack(robot);
  test_evalConstraint(robot);
  test_evalDerivatives(robot);
  testCondenseSlackAndDual(robot);
  testExpandSlackAndDual(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}