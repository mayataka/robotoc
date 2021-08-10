#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_velocity_upper_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace idocp {

class JointVelocityUpperLimitTest : public ::testing::Test {
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


void JointVelocityUpperLimitTest::testKinematics(Robot& robot) const {
  JointVelocityUpperLimit limit(robot); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::VelocityLevel);
}


void JointVelocityUpperLimitTest::testIsFeasible(Robot& robot) const {
  JointVelocityUpperLimit limit(robot); 
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  EXPECT_EQ(limit.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.v = 2*robot.jointVelocityLimit();
  EXPECT_FALSE(limit.isFeasible(robot, data, s));
}


void JointVelocityUpperLimitTest::testSetSlack(Robot& robot) const {
  JointVelocityUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc(), limit.barrier()), data_ref(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  limit.setSlack(robot, data, s);
  data_ref.slack = vmax - s.v.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityUpperLimitTest::testComputePrimalAndDualResidual(Robot& robot) const {
  JointVelocityUpperLimit limit(robot); 
  const int dimc = limit.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  data.slack.setRandom();
  data.dual.setRandom();
  auto data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  data_ref.residual = s.v.tail(dimc) - vmax + data_ref.slack;
  pdipm::ComputeComplementarySlackness(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityUpperLimitTest::testComputePrimalResidualDerivatives(Robot& robot) const {
  JointVelocityUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const auto s = SplitSolution::Random(robot);
  limit.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  limit.computePrimalResidualDerivatives(robot, data, dt, s, kkt_res);
  kkt_res_ref.lv().tail(dimc) += dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointVelocityUpperLimitTest::testCondenseSlackAndDual(Robot& robot) const {
  JointVelocityUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  limit.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  kkt_res_ref.lv().tail(dimc).array() 
      += dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qvv().diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointVelocityUpperLimitTest::testExpandSlackAndDual(Robot& robot) const {
  JointVelocityUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  limit.setSlack(robot, data, s);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  limit.expandSlackAndDual(data, s, d);
  data_ref.dslack = - d.dv().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointVelocityUpperLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlack(robot);
  testComputePrimalAndDualResidual(robot);
  testComputePrimalResidualDerivatives(robot);
  testCondenseSlackAndDual(robot);
  testExpandSlackAndDual(robot);
}

TEST_F(JointVelocityUpperLimitTest, floatingBase) {
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