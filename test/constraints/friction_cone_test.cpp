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
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
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
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.8;
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

  double barrier, dtau, mu;
  std::string fixed_base_urdf, floating_base_urdf;
};


void FrictionConeTest::testKinematics(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void FrictionConeTest::testIsFeasible(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), contact_status.max_point_contacts());
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
    for (int i=0; i<contact_status.max_point_contacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        if (FrictionCone::frictionConeResidual(robot.frictionCoefficient(i), s.f[i]) > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(limit.isFeasible(robot, data, s), feasible);
  }
}


void FrictionConeTest::testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, dtau, s);
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    data_ref.slack.coeffRef(i) = - dtau * FrictionCone::frictionConeResidual(robot.frictionCoefficient(i), s.f[i]);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.slack.isApprox(data_ref.slack));
  EXPECT_TRUE(data.dual.isApprox(data_ref.dual));
}


void FrictionConeTest::testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  KKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  KKTResidual kkt_res_ref = kkt_res;
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      Eigen::Vector3d gf;
      gf(0) = 2 * dtau * s.f[i].coeff(0);
      gf(1) = 2 * dtau * s.f[i].coeff(1);
      gf(2) = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      kkt_res_ref.lf().segment<3>(dimf_stack) += gf * data_ref.dual(i);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void FrictionConeTest::testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  ConstraintComponentData data(limit.dimc());
  data.slack.setRandom();
  data.dual.setRandom();
  ConstraintComponentData data_ref = data;
  limit.computePrimalAndDualResidual(robot, data, dtau, s);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = s.f[i].coeff(0);
      const double fy = s.f[i].coeff(1);
      const double fz = s.f[i].coeff(2);
      const double cone = fx*fx + fy*fy - mu*mu*fz*fz;
      data_ref.residual.coeffRef(i) = dtau * cone + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(data_ref.residual.isApprox(data.residual));
  EXPECT_TRUE(data_ref.duality.isApprox(data.duality));
}


void FrictionConeTest::testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, dtau, s);
  ConstraintComponentData data_ref = data;
  KKTMatrix kkt_mat(robot);
  kkt_mat.setContactStatus(contact_status);
  kkt_mat.Qff().setRandom();
  KKTResidual kkt_res(robot);
  kkt_res.setContactStatus(contact_status);
  kkt_res.lf().setRandom();
  KKTMatrix kkt_mat_ref = kkt_mat;
  KKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      const double fx = s.f[i].coeff(0);
      const double fy = s.f[i].coeff(1);
      const double fz = s.f[i].coeff(2);
      const double cone = fx*fx + fy*fy - mu*mu*fz*fz;
      data_ref.residual.coeffRef(i) = dtau * cone + data_ref.slack.coeff(i);
      data_ref.duality.coeffRef(i) = data_ref.slack.coeff(i) * data_ref.dual.coeff(i) - barrier;
      dimf_stack += 3;
    }
  }
  dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      Eigen::Vector3d gf;
      gf(0) = 2 * dtau * s.f[i].coeff(0);
      gf(1) = 2 * dtau * s.f[i].coeff(1);
      gf(2) = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      kkt_res_ref.lf().segment<3>(dimf_stack) 
          += gf * (data_ref.dual.coeff(i)*data_ref.residual.coeff(i)-data_ref.duality.coeff(i)) 
                / data_ref.slack.coeff(i);
      kkt_mat_ref.Qff().block<3, 3>(dimf_stack, dimf_stack)
          += gf * gf.transpose() * data_ref.dual.coeff(i) / data_ref.slack.coeff(i);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void FrictionConeTest::testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const {
  FrictionCone limit(robot); 
  ConstraintComponentData data(limit.dimc());
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, dtau, s);
  data.residual.setRandom();
  data.duality.setRandom();
  ConstraintComponentData data_ref = data;
  const SplitDirection d = SplitDirection::Random(robot, contact_status);
  limit.computeSlackAndDualDirection(robot, data, dtau, s, d);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      const double mu = robot.frictionCoefficient(i);
      Eigen::Vector3d gf;
      gf(0) = 2 * dtau * s.f[i].coeff(0);
      gf(1) = 2 * dtau * s.f[i].coeff(1);
      gf(2) = - 2 * dtau * mu * mu * s.f[i].coeff(2);
      data_ref.dslack.coeffRef(i) = - gf.dot(d.df().segment<3>(dimf_stack)) - data_ref.residual.coeff(i);
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


TEST_F(FrictionConeTest, frictionConeResidual) {
  const double mu = 0.8;
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
  ContactStatus contact_status(frames.size());
  contact_status.setContactStatus({false});
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  contact_status.setContactStatus({true});
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
  ContactStatus contact_status(frames.size());
  contact_status.setContactStatus({false, false, false, false});
  testKinematics(robot, contact_status);
  testIsFeasible(robot, contact_status);
  testSetSlackAndDual(robot, contact_status);
  testAugmentDualResidual(robot, contact_status);
  testComputePrimalAndDualResidual(robot, contact_status);
  testCondenseSlackAndDual(robot, contact_status);
  testComputeSlackAndDualDirection(robot, contact_status);
  std::random_device rnd;
  for (int i=0; i<contact_status.max_point_contacts(); ++i) {
    if (rnd()%2 == 0) 
      contact_status.activateContact(i);
  }
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