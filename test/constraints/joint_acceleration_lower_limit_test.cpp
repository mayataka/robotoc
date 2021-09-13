#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_acceleration_lower_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

#include "robot_factory.hpp"

namespace idocp {

class JointAccelerationLowerLimitTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const Eigen::VectorXd& amin) const;
  void testIsFeasible(Robot& robot, const Eigen::VectorXd& amin) const;
  void testSetSlack(Robot& robot, const Eigen::VectorXd& amin) const;
  void test_evalConstraint(Robot& robot, const Eigen::VectorXd& amin) const;
  void test_evalDerivatives(Robot& robot, const Eigen::VectorXd& amin) const;
  void testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const;
  void testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const;

  double barrier, dt;
};


void JointAccelerationLowerLimitTest::testKinematics(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  EXPECT_FALSE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void JointAccelerationLowerLimitTest::testIsFeasible(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  EXPECT_EQ(constr.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(constr.isFeasible(robot, data, s));
  s.a = 2*amin;
  EXPECT_FALSE(constr.isFeasible(robot, data, s));
}


void JointAccelerationLowerLimitTest::testSetSlack(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter()), data_ref(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  data_ref.slack = -amin + s.a.tail(dimc);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationLowerLimitTest::test_evalConstraint(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  auto data_ref = data;
  constr.evalConstraint(robot, data, s);
  data_ref.residual = - s.a.tail(dimc) + amin + data_ref.slack;
  pdipm::computeComplementarySlackness(barrier, data_ref);
  data_ref.log_barrier = pdipm::logBarrier(barrier, data_ref.slack);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void JointAccelerationLowerLimitTest::test_evalDerivatives(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, data, dt, s, kkt_res);
  kkt_res_ref.la.tail(dimc) -= dt * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointAccelerationLowerLimitTest::testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
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
      -= dt * (data_ref.dual.array()*data_ref.residual.array()-data_ref.cmpl.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qaa.diagonal().tail(dimc).array() 
      += dt * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointAccelerationLowerLimitTest::testExpandSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit constr(robot, amin); 
  ConstraintComponentData data(constr.dimc(), constr.barrierParameter());
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot);
  constr.setSlack(robot, data, s);
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot);
  constr.expandSlackAndDual(data, s, d);
  data_ref.dslack = d.da().tail(dimc) - data_ref.residual;
  pdipm::computeDualDirection(data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(JointAccelerationLowerLimitTest, fixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  const Eigen::VectorXd amin = Eigen::VectorXd::Constant(robot.dimv(), -10);
  testKinematics(robot, amin);
  testIsFeasible(robot, amin);
  testSetSlack(robot, amin);
  test_evalConstraint(robot, amin);
  test_evalDerivatives(robot, amin);
  testCondenseSlackAndDual(robot, amin);
  testExpandSlackAndDual(robot, amin);
}


TEST_F(JointAccelerationLowerLimitTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  const Eigen::VectorXd amin = Eigen::VectorXd::Constant(robot.dimu(), -10);
  testKinematics(robot, amin);
  testIsFeasible(robot, amin);
  testSetSlack(robot, amin);
  test_evalConstraint(robot, amin);
  test_evalDerivatives(robot, amin);
  testCondenseSlackAndDual(robot, amin);
  testExpandSlackAndDual(robot, amin);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}