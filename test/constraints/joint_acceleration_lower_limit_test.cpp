#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/joint_acceleration_lower_limit.hpp"
#include "idocp/constraints/pdipm.hpp"

namespace idocp {

class JointAccelerationLowerLimitTest : public ::testing::Test {
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

  void testKinematics(Robot& robot, const Eigen::VectorXd& amin) const;
  void testIsFeasible(Robot& robot, const Eigen::VectorXd& amin) const;
  void testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const;
  void testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amin) const;
  void testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amin) const;
  void testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const;
  void testComputeSlackAndDualDirection(Robot& robot, const Eigen::VectorXd& amin) const;

  double barrier, dtau;
  std::string fixed_base_urdf, floating_base_urdf;
};


void JointAccelerationLowerLimitTest::testKinematics(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void JointAccelerationLowerLimitTest::testIsFeasible(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), robot.dimv()-robot.dim_passive());
  SplitSolution s(robot);
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.a = 2*amin;
  EXPECT_FALSE(limit.isFeasible(robot, data, s));
}


void JointAccelerationLowerLimitTest::testSetSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  data_ref.slack = dtau * (-amin+s.a.tail(dimc));
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


void JointAccelerationLowerLimitTest::testAugmentDualResidual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  SplitKKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  kkt_res_ref.la.tail(dimc) -= dtau * data_ref.dual;
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void JointAccelerationLowerLimitTest::testComputePrimalAndDualResidual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, dtau, s);
  data_ref.residual = dtau * (- s.a.tail(dimc) + amin) + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  EXPECT_TRUE(data_ref.residual.isApprox(data.residual));
  EXPECT_TRUE(data_ref.duality.isApprox(data.duality));
}


void JointAccelerationLowerLimitTest::testCondenseSlackAndDual(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.Qaa().setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.la.setRandom();
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  data_ref.residual = dtau * (- s.a.tail(dimc) + amin) + data_ref.slack;
  pdipm::ComputeDuality(barrier, data_ref);
  kkt_res_ref.la.tail(dimc).array() 
      -= dtau * (data_ref.dual.array()*data_ref.residual.array()-data_ref.duality.array()) 
               / data_ref.slack.array();
  kkt_mat_ref.Qaa().diagonal().tail(dimc).array() 
      += dtau * dtau * data_ref.dual.array() / data_ref.slack.array();
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void JointAccelerationLowerLimitTest::testComputeSlackAndDualDirection(Robot& robot, const Eigen::VectorXd& amin) const {
  JointAccelerationLowerLimit limit(robot, amin); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot);
  limit.setSlackAndDual(robot, data, dtau, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot);
  limit.computeSlackAndDualDirection(robot, data, dtau, s, d);
  data_ref.dslack = dtau * d.da().tail(dimc) - data_ref.residual;
  pdipm::ComputeDualDirection(data_ref);
  EXPECT_TRUE(data.dslack.isApprox(data_ref.dslack));
  EXPECT_TRUE(data.ddual.isApprox(data_ref.ddual));
}


TEST_F(JointAccelerationLowerLimitTest, fixedBase) {
  Robot robot(fixed_base_urdf);
  const Eigen::VectorXd amin = Eigen::VectorXd::Constant(robot.dimv(), -10);
  testKinematics(robot, amin);
  testIsFeasible(robot, amin);
  testSetSlackAndDual(robot, amin);
  testAugmentDualResidual(robot, amin);
  testComputePrimalAndDualResidual(robot, amin);
  testCondenseSlackAndDual(robot, amin);
  testComputeSlackAndDualDirection(robot, amin);
}


TEST_F(JointAccelerationLowerLimitTest, floatingBase) {
  Robot robot(floating_base_urdf);
  const Eigen::VectorXd amin = Eigen::VectorXd::Constant(robot.dimu(), -10);
  testKinematics(robot, amin);
  testIsFeasible(robot, amin);
  testSetSlackAndDual(robot, amin);
  testAugmentDualResidual(robot, amin);
  testComputePrimalAndDualResidual(robot, amin);
  testCondenseSlackAndDual(robot, amin);
  testComputeSlackAndDualDirection(robot, amin);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}