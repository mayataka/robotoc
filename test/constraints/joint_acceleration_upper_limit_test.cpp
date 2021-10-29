#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/constraints/joint_acceleration_upper_limit.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace robotoc {

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
  void testSetSlack(Robot& robot, const Eigen::VectorXd& amax) const;
  void test_evalConstraint(Robot& robot, const Eigen::VectorXd& amax) const;
  void test_evalDerivatives(Robot& robot, const Eigen::VectorXd& amax) const;
  void testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;

  double barrier, dt;
};


void JointAccelerationUpperLimitTest::testKinematics(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  EXPECT_FALSE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void JointAccelerationUpperLimitTest::testIsFeasible(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, data, s));
  s.a = 2*amax;
  EXPECT_FALSE(constr.isFeasible(robot, data, s));
}


void JointAccelerationUpperLimitTest::testSetSlack(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter()), data_ref(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  data_ref.slack = amax - s.a.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationUpperLimitTest::test_evalConstraint(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  auto data_ref = data;
  constr.evalConstraint(robot, data, s);
  data_ref.residual = s.a.tail(dimc) - amax + data_ref.slack;
  pdipm::computeComplementarySlackness(barrier, data_ref);
  data_ref.log_barrier = pdipm::logBarrier(barrier, data_ref.slack);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationUpperLimitTest::test_evalDerivatives(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, data, dt, s, kkt_res);
  kkt_res_ref.la.tail(dimc) += dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointAccelerationUpperLimitTest::testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_mat = SplitKKTMatrix::Random(robot);
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  kkt_res_ref.la.tail(dimc).array() 
      += dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qaa.diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointAccelerationUpperLimitTest::testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit constr(robot, amax); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  constr.expandSlackAndDual(data, s, d);
  data_ref.dslack = - d.da().tail(dimc) - data_ref.residual;
  pdipm::computeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointAccelerationUpperLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimv(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlack(robot, amax);
  test_evalConstraint(robot, amax);
  test_evalDerivatives(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testExpandSlackAndDual(robot, amax);
}


TEST_F(JointAccelerationUpperLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimu(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlack(robot, amax);
  test_evalConstraint(robot, amax);
  test_evalDerivatives(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testExpandSlackAndDual(robot, amax);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}