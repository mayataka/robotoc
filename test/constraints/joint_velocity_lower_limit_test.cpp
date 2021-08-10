#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_velocity_lower_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace idocp {

class JointVelocityLowerLimitTest : public ::testing::Test {
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
  void testComputePrimalAndDualResidual(Robot& robot) const;
  void testComputePrimalResidualDerivatives(Robot& robot) const;
  void testCondenseSlackAndDual(Robot& robot) const;
  void testExpandSlackAndDual(Robot& robot) const;

  double barrier, dt;
};


void JointVelocityLowerLimitTest::testKinematics(Robot& robot) const {
  JointVelocityLowerLimit constr(robot); 
  EXPECT_FALSE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::VelocityLevel);
}


void JointVelocityLowerLimitTest::testIsFeasible(Robot& robot) const {
  JointVelocityLowerLimit constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, data, s));
  s.v = -2*robot.jointVelocityLimit();
  EXPECT_FALSE(constr.isFeasible(robot, data, s));
}


void JointVelocityLowerLimitTest::testSetSlack(Robot& robot) const {
  JointVelocityLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter()), data_ref(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmin = - robot.jointVelocityLimit();
  constr.setSlack(robot, data, s);
  data_ref.slack = -vmin + s.v.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityLowerLimitTest::testComputePrimalAndDualResidual(Robot& robot) const {
  JointVelocityLowerLimit constr(robot); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmin = - robot.jointVelocityLimit();
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  auto data_ref = data;
  constr.computePrimalAndDualResidual(robot, data, s);
  data_ref.residual = - s.v.tail(dimc) + vmin + data_ref.slack;
  pdipm::ComputeComplementarySlackness(barrier, data_ref);
  data_ref.log_barrier = pdipm::LogBarrier(barrier, data_ref.slack);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityLowerLimitTest::testComputePrimalResidualDerivatives(Robot& robot) const {
  JointVelocityLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  constr.computePrimalResidualDerivatives(robot, data, dt, s, kkt_res);
  kkt_res_ref.lv().tail(dimc) -= dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointVelocityLowerLimitTest::testCondenseSlackAndDual(Robot& robot) const {
  JointVelocityLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmin = - robot.jointVelocityLimit();
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  kkt_res_ref.lv().tail(dimc).array() 
      -= dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qvv().diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointVelocityLowerLimitTest::testExpandSlackAndDual(Robot& robot) const {
  JointVelocityLowerLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmin = - robot.jointVelocityLimit();
  constr.setSlack(robot, data, s);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  constr.expandSlackAndDual(data, s, d);
  data_ref.dslack = d.dv().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointVelocityLowerLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlack(robot);
  testComputePrimalAndDualResidual(robot);
  testComputePrimalResidualDerivatives(robot);
  testCondenseSlackAndDual(robot);
  testExpandSlackAndDual(robot);
}


TEST_F(JointVelocityLowerLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlack(robot);
  testComputePrimalAndDualResidual(robot);
  testComputePrimalResidualDerivatives(robot);
  testCondenseSlackAndDual(robot);
  testExpandSlackAndDual(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}