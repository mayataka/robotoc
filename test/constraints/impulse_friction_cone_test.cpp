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
#include "idocp/constraints/impulse_friction_cone.hpp"

namespace idocp {

class ImpulseFrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
    mu = 0.7;
    fraction_to_boundary_rate = 0.995;
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

  double barrier, mu, fraction_to_boundary_rate;
  std::string fixed_base_urdf, floating_base_urdf;
};


void ImpulseFrictionConeTest::testKinematics(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void ImpulseFrictionConeTest::testIsFeasible(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), 2*impulse_status.maxPointContacts());
  ImpulseSplitSolution s(robot);
  s.setImpulseStatus(impulse_status);
  s.f_stack().setZero();
  s.set_f_vector();
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  if (impulse_status.hasActiveImpulse()) {
    s.f_stack().setConstant(1.0);
    s.set_f_vector();
    EXPECT_FALSE(limit.isFeasible(robot, data, s));
    s.f_stack().setRandom();
    s.set_f_vector();
    bool feasible = true;
    for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
      if (impulse_status.isImpulseActive(i)) {
        if (ImpulseFrictionCone::normalForceResidual(s.f[i]) > 0) {
          feasible = false;
        }
        if (ImpulseFrictionCone::frictionConeResidual(mu, s.f[i]) > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(limit.isFeasible(robot, data, s), feasible);
  }
}


void ImpulseFrictionConeTest::testSetSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  limit.allocateExtraData(data);
  limit.allocateExtraData(data_ref);
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    data_ref.slack(2*i)   = - ImpulseFrictionCone::normalForceResidual(s.f[i]);
    data_ref.slack(2*i+1) = - ImpulseFrictionCone::frictionConeResidual(mu, s.f[i]);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void ImpulseFrictionConeTest::testAugmentDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
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
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().coeffRef(dimf_stack+2) -= data_ref.dual(2*i);
      Eigen::Vector3d gf;
      gf(0) = 2 * s.f[i].coeff(0);
      gf(1) = 2 * s.f[i].coeff(1);
      gf(2) = - 2 * mu * mu * s.f[i].coeff(2);
      EXPECT_TRUE(gf.isApprox(data.r[i]));
      kkt_res_ref.lf().segment<3>(dimf_stack) += gf * data_ref.dual(2*i+1);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void ImpulseFrictionConeTest::testComputePrimalAndDualResidual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.residual(2*i)   = ImpulseFrictionCone::normalForceResidual(s.f[i]) + data_ref.slack(2*i);
      data_ref.residual(2*i+1) = ImpulseFrictionCone::frictionConeResidual(mu, s.f[i]) + data_ref.slack(2*i+1);
      data_ref.duality(2*i)   = data_ref.slack(2*i)   * data_ref.dual(2*i)   - barrier;
      data_ref.duality(2*i+1) = data_ref.slack(2*i+1) * data_ref.dual(2*i+1) - barrier;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void ImpulseFrictionConeTest::testCondenseSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
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
  limit.augmentDualResidual(robot, data, s, kkt_res);
  ImpulseSplitKKTMatrix kkt_mat_ref = kkt_mat;
  ImpulseSplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, s, kkt_mat, kkt_res);
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.residual(2*i)   = ImpulseFrictionCone::normalForceResidual(s.f[i]) + data_ref.slack(2*i);
      data_ref.residual(2*i+1) = ImpulseFrictionCone::frictionConeResidual(mu, s.f[i]) + data_ref.slack(2*i+1);
      data_ref.duality(2*i)   = data_ref.slack(2*i)   * data_ref.dual(2*i)   - barrier;
      data_ref.duality(2*i+1) = data_ref.slack(2*i+1) * data_ref.dual(2*i+1) - barrier;
    }
  }
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      kkt_res_ref.lf().coeffRef(dimf_stack+2) 
          -= (data_ref.dual(2*i)*data_ref.residual(2*i)-data_ref.duality(2*i)) 
                  / data_ref.slack(2*i);
      kkt_mat_ref.Qff().coeffRef(dimf_stack+2, dimf_stack+2)
          += data_ref.dual(2*i) / data_ref.slack(2*i);
      Eigen::Vector3d gf;
      gf(0) = 2 * s.f[i].coeff(0);
      gf(1) = 2 * s.f[i].coeff(1);
      gf(2) = - 2 * mu * mu * s.f[i].coeff(2);
      EXPECT_TRUE(gf.isApprox(data.r[i]));
      kkt_res_ref.lf().segment<3>(dimf_stack) 
          += gf * (data_ref.dual(2*i+1)*data_ref.residual(2*i+1)-data_ref.duality(2*i+1)) 
                  / data_ref.slack(2*i+1);
      kkt_mat_ref.Qff().block<3, 3>(dimf_stack, dimf_stack)
          += gf * gf.transpose() * data_ref.dual(2*i+1) / data_ref.slack(2*i+1);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void ImpulseFrictionConeTest::testComputeSlackAndDualDirection(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  const int dimc = limit.dimc();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot, impulse_status);
  limit.setSlackAndDual(robot, data, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ImpulseSplitKKTMatrix kkt_mat(robot);
  kkt_mat.setImpulseStatus(impulse_status);
  kkt_mat.Qff().setRandom();
  ImpulseSplitKKTResidual kkt_res(robot);
  kkt_res.setImpulseStatus(impulse_status);
  kkt_res.lf().setRandom();
  limit.augmentDualResidual(robot, data, s, kkt_res);
  limit.condenseSlackAndDual(robot, data, s, kkt_mat, kkt_res);
  ConstraintComponentData data_ref = data;
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot, impulse_status);
  limit.computeSlackAndDualDirection(robot, data, s, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<impulse_status.maxPointContacts(); ++i) {
    if (impulse_status.isImpulseActive(i)) {
      data_ref.dslack(2*i) = d.df().coeff(dimf_stack+2) - data_ref.residual(2*i);
      data_ref.ddual(2*i)  = - (data_ref.dual(2*i)*data_ref.dslack(2*i)+data_ref.duality(2*i))
                                      / data_ref.slack(2*i);
      Eigen::Vector3d gf;
      gf(0) = 2 * s.f[i].coeff(0);
      gf(1) = 2 * s.f[i].coeff(1);
      gf(2) = - 2 * mu * mu * s.f[i].coeff(2);
      EXPECT_TRUE(gf.isApprox(data.r[i]));
      data_ref.dslack(2*i+1) = - gf.dot(d.df().segment<3>(dimf_stack)) - data_ref.residual(2*i+1);
      data_ref.ddual(2*i+1) = - (data_ref.dual(2*i+1)*data_ref.dslack(2*i+1)+data_ref.duality(2*i+1))
                                      / data_ref.slack(2*i+1);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(ImpulseFrictionConeTest, normalForceResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  const double fz = f(2);
  const double normal_ref = - fz;
  EXPECT_DOUBLE_EQ(normal_ref, ImpulseFrictionCone::normalForceResidual(f));
}


TEST_F(ImpulseFrictionConeTest, frictionConeResidual) {
  const double mu = 0.8;
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  const double fx = f(0);
  const double fy = f(1);
  const double fz = f(2);
  const double cone_ref = fx*fx + fy*fy - mu*mu*fz*fz;
  EXPECT_DOUBLE_EQ(cone_ref, ImpulseFrictionCone::frictionConeResidual(mu, f));
}


TEST_F(ImpulseFrictionConeTest, fixedBase) {
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


TEST_F(ImpulseFrictionConeTest, floatingBase) {
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