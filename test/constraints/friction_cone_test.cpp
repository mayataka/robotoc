#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/constraints/pdipm.hpp"
#include "idocp/constraints/friction_cone.hpp"

namespace idocp {

class FrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.7;
    fraction_to_boundary_rate = 0.995;
  }

  virtual void TearDown() {
  }

  void testKinematics(Robot& robot, const ContactStatus& contact_status) const;
  void testIsFeasible(Robot& robot, const ContactStatus& contact_status) const;
  void testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;
  void testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const;
  void testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const;
  void testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;
  void testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const;

  double barrier, dt, mu, fraction_to_boundary_rate;
  std::string fixed_base_urdf, floating_base_urdf;
};


void FrictionConeTest::testKinematics(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void FrictionConeTest::testIsFeasible(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), 2*contact_status.maxPointContacts());
  SplitSolution s(robot);
  s.setContactStatus(contact_status);
  s.f_stack().setZero();
  s.set_f_vector();
  EXPECT_TRUE(limit.isFeasible(robot, data, s));
  if (contact_status.hasActiveContacts()) {
    s.f_stack().setConstant(1.0);
    s.set_f_vector();
    EXPECT_FALSE(limit.isFeasible(robot, data, s));
    s.f_stack().setRandom();
    s.set_f_vector();
    bool feasible = true;
    for (int i=0; i<contact_status.maxPointContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        if (FrictionCone::normalForceResidual(s.f[i]) > 0) {
          feasible = false;
        }
        if (FrictionCone::frictionConeResidual(mu, s.f[i]) > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(limit.isFeasible(robot, data, s), feasible);
  }
}


void FrictionConeTest::testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  limit.allocateExtraData(data);
  limit.allocateExtraData(data_ref);
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    data_ref.slack(2*i)   = - FrictionCone::normalForceResidual(s.f[i]);
    data_ref.slack(2*i+1) = - FrictionCone::frictionConeResidual(mu, s.f[i]);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_res_ref.lf().coeffRef(dimf_stack+2) -= dt * data_ref.dual(2*i);
      Eigen::Vector3d gf;
      gf(0) = 2 * s.f[i].coeff(0);
      gf(1) = 2 * s.f[i].coeff(1);
      gf(2) = - 2 * mu * mu * s.f[i].coeff(2);
      EXPECT_TRUE(gf.isApprox(data.r[i]));
      kkt_res_ref.lf().segment<3>(dimf_stack) += dt * gf * data_ref.dual(2*i+1);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void FrictionConeTest::testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, s);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      data_ref.residual(2*i)   = FrictionCone::normalForceResidual(s.f[i]) + data_ref.slack(2*i);
      data_ref.residual(2*i+1) = FrictionCone::frictionConeResidual(mu, s.f[i]) + data_ref.slack(2*i+1);
      data_ref.duality(2*i)   = data_ref.slack(2*i)   * data_ref.dual(2*i)   - barrier;
      data_ref.duality(2*i+1) = data_ref.slack(2*i+1) * data_ref.dual(2*i+1) - barrier;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void FrictionConeTest::testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  ConstraintComponentData data_ref = data;
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  kkt_mat.Qff().setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      data_ref.residual(2*i)   = FrictionCone::normalForceResidual(s.f[i]) + data_ref.slack(2*i);
      data_ref.residual(2*i+1) = FrictionCone::frictionConeResidual(mu, s.f[i]) + data_ref.slack(2*i+1);
      data_ref.duality(2*i)   = data_ref.slack(2*i)   * data_ref.dual(2*i)   - barrier;
      data_ref.duality(2*i+1) = data_ref.slack(2*i+1) * data_ref.dual(2*i+1) - barrier;
    }
  }
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_res_ref.lf().coeffRef(dimf_stack+2) 
          -= dt * (data_ref.dual(2*i)*data_ref.residual(2*i)-data_ref.duality(2*i)) 
                  / data_ref.slack(2*i);
      kkt_mat_ref.Qff().coeffRef(dimf_stack+2, dimf_stack+2)
          += dt * data_ref.dual(2*i) / data_ref.slack(2*i);
      Eigen::Vector3d gf;
      gf(0) = 2 * s.f[i].coeff(0);
      gf(1) = 2 * s.f[i].coeff(1);
      gf(2) = - 2 * mu * mu * s.f[i].coeff(2);
      EXPECT_TRUE(gf.isApprox(data.r[i]));
      kkt_res_ref.lf().segment<3>(dimf_stack) 
          += dt * gf * (data_ref.dual(2*i+1)*data_ref.residual(2*i+1)-data_ref.duality(2*i+1)) 
                  / data_ref.slack(2*i+1);
      kkt_mat_ref.Qff().block<3, 3>(dimf_stack, dimf_stack)
          += dt * gf * gf.transpose() * data_ref.dual(2*i+1) / data_ref.slack(2*i+1);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void FrictionConeTest::testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  limit.allocateExtraData(data);
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  data.residual.setRandom();
  data.duality.setRandom();
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
  SplitKKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  kkt_mat.Qff().setRandom();
  SplitKKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  limit.augmentDualResidual(robot, data, dt, s, kkt_res);
  limit.condenseSlackAndDual(robot, data, dt, s, kkt_mat, kkt_res);
  ConstraintComponentData data_ref = data;
  limit.computeSlackAndDualDirection(robot, data, s, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
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


TEST_F(FrictionConeTest, normalForceResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  const double fz = f(2);
  const double normal_ref = - fz;
  EXPECT_DOUBLE_EQ(normal_ref, FrictionCone::normalForceResidual(f));
}


TEST_F(FrictionConeTest, frictionConeResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  const double fx = f(0);
  const double fy = f(1);
  const double fz = f(2);
  const double cone_ref = fx*fx + fy*fy - mu*mu*fz*fz;
  EXPECT_DOUBLE_EQ(cone_ref, FrictionCone::frictionConeResidual(mu, f));
}


TEST_F(FrictionConeTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  contact_status.activateContact(0);
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
}


TEST_F(FrictionConeTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  auto contact_status = robot.createContactStatus();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  contact_status.setRandom();
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}