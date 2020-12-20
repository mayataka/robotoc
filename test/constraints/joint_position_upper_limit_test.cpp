#include <string>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_position_upper_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

namespace idocp {

class JointPositionUpperLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot) const;
  void testIsFeasible(Robot& robot) const;
  void testSetSlackAndDual(Robot& robot) const;
  void testAugmentDualResidual(Robot& robot) const;
  void testComputePrimalAndDualResidual(Robot& robot) const;
  void testCondenseSlackAndDual(Robot& robot) const;
  void testComputeSlackAndDualDirection(Robot& robot) const;

  double barrier, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


void JointPositionUpperLimitTest::testKinematics(Robot& robot) const {
  JointPositionUpperLimit limit(robot); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::PositionLevel);
}


void JointPositionUpperLimitTest::testIsFeasible(Robot& robot) const {
  JointPositionUpperLimit limit(robot); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.q = 2*robot.upperJointPositionLimit();
  EXPECT_FALSE(limit.isFeasible(robot, data, s));
}


void JointPositionUpperLimitTest::testSetSlackAndDual(Robot& robot) const {
  JointPositionUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  limit.setSlackAndDual(robot, data, s);
  data_ref.slack = qmax - s.q.tail(dimc);
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


void JointPositionUpperLimitTest::testAugmentDualResidual(Robot& robot) const {
  JointPositionUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTResidual kkt_res(robot);
  kkt_res.lq().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  kkt_res_ref.lq().tail(dimc) += dtau * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointPositionUpperLimitTest::testComputePrimalAndDualResidual(Robot& robot) const {
  JointPositionUpperLimit limit(robot); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  data_ref.residual = s.q.tail(dimc) - qmax + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  EXPECT_TRUE(data_ref.residual.isApprox(data.residual));
  EXPECT_TRUE(data_ref.duality.isApprox(data.duality));
}


void JointPositionUpperLimitTest::testCondenseSlackAndDual(Robot& robot) const {
  JointPositionUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.Qqq().setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.lq().setRandom();
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  data_ref.residual = s.q.tail(dimc) - qmax + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  kkt_res_ref.lq().tail(dimc).array() 
      += dtau * (data_ref.dual.array()*data_ref.residual.array()-data_ref.duality.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qqq().diagonal().tail(dimc).array() 
      += dtau * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointPositionUpperLimitTest::testComputeSlackAndDualDirection(Robot& robot) const {
  JointPositionUpperLimit limit(robot);
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  const Eigen::VectorXd qmax = robot.upperJointPositionLimit();
  limit.setSlackAndDual(robot, data, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot);
  limit.computeSlackAndDualDirection(robot, data, s, d);
  data_ref.dslack = - d.dq().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.dslack.isApprox(data_ref.dslack));
  EXPECT_TRUE(data.ddual.isApprox(data_ref.ddual));
}


TEST_F(JointPositionUpperLimitTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlackAndDual(robot);
  testAugmentDualResidual(robot);
  testComputePrimalAndDualResidual(robot);
  testCondenseSlackAndDual(robot);
  testComputeSlackAndDualDirection(robot);
}


TEST_F(JointPositionUpperLimitTest, floatingBase) {
  Robot robot(floating_base_urdf);
  testKinematics(robot);
  testIsFeasible(robot);
  testSetSlackAndDual(robot);
  testAugmentDualResidual(robot);
  testComputePrimalAndDualResidual(robot);
  testCondenseSlackAndDual(robot);
  testComputeSlackAndDualDirection(robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}