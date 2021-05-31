#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_acceleration_upper_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace idocp {

class JointAccelerationUpperLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const Eigen::VectorXd& amax) const;
  void testIsFeasible(Robot& robot, const Eigen::VectorXd& amax) const;
  void testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;

  double barrier, dt;
};


void JointAccelerationUpperLimitTest::testKinematics(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void JointAccelerationUpperLimitTest::testIsFeasible(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  EXPECT_EQ(limit.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.a = 2*amax;
  EXPECT_FALSE(limit.isFeasible(robot, data, s));
}


void JointAccelerationUpperLimitTest::testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc(), limit.barrier()), data_ref(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, s);
  data_ref.slack = amax - s.a.tail(dimc);
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationUpperLimitTest::testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  kkt_res_ref.la.tail(dimc) += dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointAccelerationUpperLimitTest::testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  data_ref.residual = s.a.tail(dimc) - amax + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationUpperLimitTest::testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.Qaa.setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  data_ref.residual = s.a.tail(dimc) - amax + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  kkt_res_ref.la.tail(dimc).array() 
      += dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.duality.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qaa.diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointAccelerationUpperLimitTest::testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc(), limit.barrier());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot);
  limit.expandSlackAndDual(data, s, d);
  data_ref.dslack = - d.da().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointAccelerationUpperLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimv(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlackAndDual(robot, amax);
  testAugmentDualResidual(robot, amax);
  testComputePrimalAndDualResidual(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testExpandSlackAndDual(robot, amax);
}


TEST_F(JointAccelerationUpperLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimu(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlackAndDual(robot, amax);
  testAugmentDualResidual(robot, amax);
  testComputePrimalAndDualResidual(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testExpandSlackAndDual(robot, amax);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}