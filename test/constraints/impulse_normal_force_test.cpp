#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/impulse_normal_force.hpp"

namespace idocp {

class ImpulseNormalForceTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testIsFeasible(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testSetSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testAugmentDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testComputePrimalAndDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testCondenseSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const;
  void testComputeSlackAndDualDirection(Robot& robot, const ImpulseStatus& impulse_status) const;

  double barrier;
  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseNormalForceTest::testKinematics(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void ImpulseNormalForceTest::testIsFeasible(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), impulse_status.max_point_contacts());
  ImpulseSplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  s.f_stack().setZero();
  s.set_f_vector();
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  s.f_stack().setConstant(1.0);
  s.set_f_vector();
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  if (impulse_status.hasActiveImpulse()) {
    s.f_stack().setConstant(-1.0);
    s.set_f_vector();
    EXPECT_FALSE(limit.isFeasible(robot, data, s));
  }
}


void ImpulseNormalForceTest::testSetSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    data_ref.slack.coeffRef(i) = s.f[i].coeff(2);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


void ImpulseNormalForceTest::testAugmentDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_res.lf().setRandom();
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, s, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).coeffRef(2) -= data_ref.dual.coeff(i);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void ImpulseNormalForceTest::testComputePrimalAndDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.residual.coeffRef(i) = - s.f[i].coeff(2) + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
    }
  }
  EXPECT_TRUE(data_ref.residual.isApprox(data.residual));
  EXPECT_TRUE(data_ref.duality.isApprox(data.duality));
}


void ImpulseNormalForceTest::testCondenseSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  ImpulseSplitKKTMatrix kkt_mat(robot);
  kkt_mat.setImpulseStatus(impulse_status);
  kkt_mat.Qff().setRandom();
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_res.lf().setRandom();
  ImpulseSplitKKTMatrix kkt_mat_ref = kkt_mat;
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, s, kkt_mat, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.residual.coeffRef(i) = - s.f[i].coeff(2) + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack).coeffRef(2) 
          -= (data_ref.dual.coeff(i)*data_ref.residual.coeff(i)-data_ref.duality.coeff(i)) 
             / data_ref.slack.coeff(i);
      kkt_mat_ref.Qff().diagonal().segment<3>(dimf_stack).coeffRef(2)
          += data_ref.dual.coeff(i) / data_ref.slack.coeff(i);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void ImpulseNormalForceTest::testComputeSlackAndDualDirection(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseNormalForce limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  limit.computeSlackAndDualDirection(robot, data, s, d);
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.max_point_contacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.dslack.coeffRef(i) = d.df().segment<3>(dimf_stack).coeff(2) - data_ref.residual.coeff(i);
      data_ref.ddual.coeffRef(i) = - (data_ref.dual.coeff(i)*data_ref.dslack.coeff(i)+data_ref.duality.coeff(i))
                                      / data_ref.slack.coeff(i);
      dimf_stack += 3;
    }
    else {
      data_ref.slack.coeffRef(i) = 1.0;
      data_ref.dslack.coeffRef(i) = 1.0;
      data_ref.dual.coeffRef(i) = 1.0;
      data_ref.ddual.coeffRef(i) = 1.0;
    }
  }
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
  EXPECT_TRUE(data.dslack.isApprox(data_ref.dslack));
  EXPECT_TRUE(data.ddual.isApprox(data_ref.ddual));
}


TEST_F(ImpulseNormalForceTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  impulse_status.setImpulseStatus({false});
  testKinematics(robot, impulse_status);
  testIsFeasible(robot, impulse_status);
  testSetSlackAndDual(robot, impulse_status);
  testAugmentDualResidual(robot, impulse_status);
  testComputePrimalAndDualResidual(robot, impulse_status);
  testCondenseSlackAndDual(robot, impulse_status);
  testComputeSlackAndDualDirection(robot, impulse_status);
  impulse_status.setImpulseStatus({true});
  testKinematics(robot, impulse_status);
  testIsFeasible(robot, impulse_status);
  testSetSlackAndDual(robot, impulse_status);
  testAugmentDualResidual(robot, impulse_status);
  testComputePrimalAndDualResidual(robot, impulse_status);
  testCondenseSlackAndDual(robot, impulse_status);
  testComputeSlackAndDualDirection(robot, impulse_status);
}


TEST_F(ImpulseNormalForceTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  ImpulseStatus impulse_status = robot.createImpulseStatus();
  impulse_status.setImpulseStatus({false, false, false, false});
  testKinematics(robot, impulse_status);
  testIsFeasible(robot, impulse_status);
  testSetSlackAndDual(robot, impulse_status);
  testAugmentDualResidual(robot, impulse_status);
  testComputePrimalAndDualResidual(robot, impulse_status);
  testCondenseSlackAndDual(robot, impulse_status);
  testComputeSlackAndDualDirection(robot, impulse_status);
  impulse_status.setRandom();
  testKinematics(robot, impulse_status);
  testIsFeasible(robot, impulse_status);
  testSetSlackAndDual(robot, impulse_status);
  testAugmentDualResidual(robot, impulse_status);
  testComputePrimalAndDualResidual(robot, impulse_status);
  testCondenseSlackAndDual(robot, impulse_status);
  testComputeSlackAndDualDirection(robot, impulse_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}