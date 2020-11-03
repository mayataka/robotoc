#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/constraints/joint_acceleration_upper_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

namespace idocp {

class JointAccelerationUpperLimitTest : public ::testing::Test {
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

  void testKinematics(Robot& robot, const Eigen::VectorXd& amax) const;
  void testIsFeasible(Robot& robot, const Eigen::VectorXd& amax) const;
  void testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const;
  void testComputeSlackAndDualDirection(Robot& robot, const Eigen::VectorXd& amax) const;

  double barrier, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


void JointAccelerationUpperLimitTest::testKinematics(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void JointAccelerationUpperLimitTest::testIsFeasible(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.a = 2*amax;
  EXPECT_FALSE(limit.isFeasible(robot, data, s));
}


void JointAccelerationUpperLimitTest::testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  data_ref.slack = dtau * (amax-s.a.tail(dimc));
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


void JointAccelerationUpperLimitTest::testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  KKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  KKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  kkt_res_ref.la.tail(dimc) += dtau * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointAccelerationUpperLimitTest::testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, dtau, s);
  data_ref.residual = dtau * (s.a.tail(dimc) - amax) + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  EXPECT_TRUE(data_ref.residual.isApprox(data.residual));
  EXPECT_TRUE(data_ref.duality.isApprox(data.duality));
}


void JointAccelerationUpperLimitTest::testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  KKTMatrix kkt_mat(robot);
  kkt_mat.Qaa().setRandom();
  KKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  KKTMatrix kkt_mat_ref = kkt_mat;
  KKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  data_ref.residual = dtau * (s.a.tail(dimc)-amax) + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  kkt_res_ref.la.tail(dimc).array() 
      += dtau * (data_ref.dual.array()*data_ref.residual.array()-data_ref.duality.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qaa().diagonal().tail(dimc).array() 
      += dtau * dtau * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointAccelerationUpperLimitTest::testComputeSlackAndDualDirection(Robot& robot, const Eigen::VectorXd& amax) const {
  JointAccelerationUpperLimit limit(robot, amax); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot);
  limit.computeSlackAndDualDirection(robot, data, dtau, s, d);
  data_ref.dslack = - dtau * d.da().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.dslack.isApprox(data_ref.dslack));
  EXPECT_TRUE(data.ddual.isApprox(data_ref.ddual));
}


TEST_F(JointAccelerationUpperLimitTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimv(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlackAndDual(robot, amax);
  testAugmentDualResidual(robot, amax);
  testComputePrimalAndDualResidual(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testComputeSlackAndDualDirection(robot, amax);
}


TEST_F(JointAccelerationUpperLimitTest, floatingBase) {
  Robot robot(floating_base_urdf);
  const Eigen::VectorXd amax = Eigen::VectorXd::Constant(robot.dimu(), 10);
  testKinematics(robot, amax);
  testIsFeasible(robot, amax);
  testSetSlackAndDual(robot, amax);
  testAugmentDualResidual(robot, amax);
  testComputePrimalAndDualResidual(robot, amax);
  testCondenseSlackAndDual(robot, amax);
  testComputeSlackAndDualDirection(robot, amax);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}