#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class JointVelocityUpperLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier_param = 1.0e-03;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    grid_info = GridInfo::Random();
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

  double barrier_param, dt;
  GridInfo grid_info;
};


void JointVelocityUpperLimitTest::test_kinematics(Robot& robot) const {
  JointVelocityUpperLimit constr(robot); 
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::VelocityLevel);
}


void JointVelocityUpperLimitTest::test_isFeasible(Robot& robot) const {
  JointVelocityUpperLimit constr(robot); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  const auto contact_status = robot.createContactStatus();
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, contact_status, grid_info, s, data));
  s.v = 2*robot.jointVelocityLimit();
  EXPECT_FALSE(constr.isFeasible(robot, contact_status, grid_info, s, data));
}


void JointVelocityUpperLimitTest::test_setSlack(Robot& robot) const {
  JointVelocityUpperLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam()), data_ref(constr.dimc(), constr.getBarrierParam());
  const int dimc = constr.dimc();
  const auto contact_status = robot.createContactStatus();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  constr.setSlack(robot, contact_status, grid_info, s, data);
  data_ref.slack = vmax - s.v.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityUpperLimitTest::test_evalConstraint(Robot& robot) const {
  JointVelocityUpperLimit constr(robot); 
  const int dimc = constr.dimc();
  const auto contact_status = robot.createContactStatus();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  auto data_ref = data;
  constr.evalConstraint(robot, contact_status, grid_info, s, data);
  data_ref.residual = s.v.tail(dimc) - vmax + data_ref.slack;
  pdipm::computeComplementarySlackness(barrier_param, data_ref);
  data_ref.log_barrier = pdipm::logBarrier(barrier_param, data_ref.slack);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointVelocityUpperLimitTest::test_evalDerivatives(Robot& robot) const {
  JointVelocityUpperLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  const int dimc = constr.dimc();
  const auto contact_status = robot.createContactStatus();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, contact_status, grid_info, s, data);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, contact_status, grid_info, s, data, kkt_res);
  kkt_res_ref.lv().tail(dimc) += data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointVelocityUpperLimitTest::test_condenseSlackAndDual(Robot& robot) const {
  JointVelocityUpperLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  const int dimc = constr.dimc();
  const auto contact_status = robot.createContactStatus();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  constr.setSlack(robot, contact_status, grid_info, s, data);
  auto data_ref = data;
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(contact_status, grid_info, data, kkt_mat, kkt_res);
  kkt_res_ref.lv().tail(dimc).array() 
      += (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qvv().diagonal().tail(dimc).array() 
      += data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointVelocityUpperLimitTest::test_expandSlackAndDual(Robot& robot) const {
  JointVelocityUpperLimit constr(robot);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  const int dimc = constr.dimc();
  const auto contact_status = robot.createContactStatus();
  const auto s = SplitSolution::Random(robot);
  const Eigen::VectorXd vmax = robot.jointVelocityLimit();
  constr.setSlack(robot, contact_status, grid_info, s, data);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  constr.expandSlackAndDual(contact_status, grid_info, d, data);
  data_ref.dslack = - d.dv().tail(dimc) - data_ref.residual;
  pdipm::computeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointVelocityUpperLimitTest, fixedBase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  test_kinematics(robot);
  test_isFeasible(robot);
  test_setSlack(robot);
  test_evalConstraint(robot);
  test_evalDerivatives(robot);
  test_condenseSlackAndDual(robot);
  test_expandSlackAndDual(robot);
}

TEST_F(JointVelocityUpperLimitTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
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