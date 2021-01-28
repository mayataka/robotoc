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
#include "idocp/constraints/linearized_friction_cone.hpp"

namespace idocp {

class LinearizedFrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    barrier = 1.0e-04;
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.7;
    fraction_to_boundary_rate = 0.995;
    Jac = Eigen::MatrixXd::Zero(5, 3);
    Jac <<  0,  0, -1,
            1,  0, -mu/std::sqrt(2),
           -1,  0, -mu/std::sqrt(2),
            0,  1, -mu/std::sqrt(2),
            0, -1, -mu/std::sqrt(2);
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

  double barrier, dtau, mu, fraction_to_boundary_rate;
  std::string fixed_base_urdf, floating_base_urdf;
  Eigen::MatrixXd Jac;
};


void LinearizedFrictionConeTest::testKinematics(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
  EXPECT_FALSE(limit.useKinematics());
  EXPECT_TRUE(limit.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void LinearizedFrictionConeTest::testIsFeasible(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc());
  EXPECT_EQ(limit.dimc(), 5*contact_status.maxPointContacts());
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
        Eigen::VectorXd res(5);
        LinearizedFrictionCone::frictionConeResidual(mu, s.f[i], res);
        if (res.maxCoeff() > 0) {
          feasible = false;
        }
      }
    }
    EXPECT_EQ(limit.isFeasible(robot, data, s), feasible);
  }
}


void LinearizedFrictionConeTest::testSetSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
  ConstraintComponentData data(limit.dimc()), data_ref(limit.dimc());
  limit.allocateExtraData(data);
  limit.allocateExtraData(data_ref);
  const int dimc = limit.dimc();
  const SplitSolution s = SplitSolution::Random(robot, contact_status);
  limit.setSlackAndDual(robot, data, s);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    LinearizedFrictionCone::frictionConeResidual(mu, s.f[i], data_ref.residual.segment(5*i, 5));
    data_ref.slack.segment(5*i, 5) = - data_ref.residual.segment(5*i, 5);
  }
  pdipm::SetSlackAndDualPositive(barrier, data_ref);
  EXPECT_TRUE(data.isApprox(data_ref));
}


void LinearizedFrictionConeTest::testAugmentDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
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
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      kkt_res_ref.lf().segment<3>(dimf_stack) += dtau * Jac.transpose() * data_ref.dual.segment(5*i, 5);
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void LinearizedFrictionConeTest::testComputePrimalAndDualResidual(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
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
      LinearizedFrictionCone::frictionConeResidual(mu, s.f[i], data_ref.residual.segment(5*i, 5));
      for (int j=0; j<5; ++j) {
        data_ref.duality(5*i+j) = data_ref.slack(5*i+j) * data_ref.dual(5*i+j) - barrier;
      }
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void LinearizedFrictionConeTest::testCondenseSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
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
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  SplitKKTMatrix kkt_mat_ref = kkt_mat;
  SplitKKTResidual kkt_res_ref = kkt_res;
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      LinearizedFrictionCone::frictionConeResidual(mu, s.f[i], data_ref.residual.segment(5*i, 5));
      for (int j=0; j<5; ++j) {
        data_ref.duality(5*i+j) = data_ref.slack(5*i+j) * data_ref.dual(5*i+j) - barrier;
      }
    }
  }
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      Eigen::VectorXd vec(5);
      vec.array() = (data_ref.dual.segment(5*i, 5).array()*data_ref.residual.segment(5*i, 5).array()-data_ref.duality.segment(5*i, 5).array()) 
                  / data_ref.slack.segment(5*i, 5).array();
      kkt_res_ref.lf().segment(dimf_stack, 3)
          += dtau * Jac.transpose() * vec;
      vec.array() = data_ref.dual.segment(5*i, 5).array() / data_ref.slack.segment(5*i, 5).array();
      kkt_mat_ref.Qff().block(dimf_stack, dimf_stack, 3, 3)
          += dtau * Jac.transpose() * vec.asDiagonal() * Jac;
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void LinearizedFrictionConeTest::testComputeSlackAndDualDirection(Robot& robot, const ContactStatus& contact_status) const {
  LinearizedFrictionCone limit(robot, mu); 
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
  limit.augmentDualResidual(robot, data, dtau, s, kkt_res);
  limit.condenseSlackAndDual(robot, data, dtau, s, kkt_mat, kkt_res);
  ConstraintComponentData data_ref = data;
  limit.computeSlackAndDualDirection(robot, data, s, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  for (int i=0; i<contact_status.maxPointContacts(); ++i) {
    if (contact_status.isContactActive(i)) {
      data_ref.dslack.segment(5*i, 5) = - Jac * d.df().segment<3>(dimf_stack) - data_ref.residual.segment(5*i, 5);
      for (int j=0; j<5; ++j) {
        data_ref.ddual(5*i+j) = - (data_ref.dual(5*i+j)*data_ref.dslack(5*i+j)+data_ref.duality(5*i+j))
                                        / data_ref.slack(5*i+j);
      }
      dimf_stack += 3;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(LinearizedFrictionConeTest, frictionConeResidual) {
  const Eigen::Vector3d f = Eigen::Vector3d::Random();
  Eigen::VectorXd res(5);
  LinearizedFrictionCone::frictionConeResidual(mu, f, res);
  Eigen::VectorXd res_ref(5);
  const double fx = f(0);
  const double fy = f(1);
  const double fz = f(2);
  res_ref(0) = - fz;
  res_ref(1) =   fx - mu * fz / std::sqrt(2);
  res_ref(2) = - fx - mu * fz / std::sqrt(2);
  res_ref(3) =   fy - mu * fz / std::sqrt(2);
  res_ref(4) = - fy - mu * fz / std::sqrt(2);
  EXPECT_TRUE(res.isApprox(res_ref));
}


TEST_F(LinearizedFrictionConeTest, fixedBase) {
  const std::vector<int> frames = {18};
  Robot robot(fixed_base_urdf, frames);
  ContactStatus contact_status = robot.createContactStatus();
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


TEST_F(LinearizedFrictionConeTest, floatingBase) {
  const std::vector<int> frames = {14, 24, 34, 44};
  Robot robot(floating_base_urdf, frames);
  ContactStatus contact_status = robot.createContactStatus();
  contact_status.setContactStatus({false, false, false, false});
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