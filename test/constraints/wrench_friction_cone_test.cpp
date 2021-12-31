#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/constraints/pdipm.hpp"
#include "robotoc/constraints/wrench_friction_cone.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class WrenchFrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier = 1.0e-03;
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
    mu = 0.7;
    X = 0.1;
    Y = 0.05;
    fraction_to_boundary_rule = 0.995;
    cone.resize(17, 6);
    cone.setZero();
    cone <<  0,  0, -1, 0, 0, 0,
            -1,  0, -mu, 0, 0, 0,
              1,  0, -mu, 0, 0, 0,
              0, -1, -mu, 0, 0, 0,
              0,  1, -mu, 0, 0, 0,
              0,  0, -Y, -1,  0, 0,
              0,  0, -Y,  1,  0, 0,
              0,  0, -X,  0, -1, 0,
              0,  0, -X,  0,  1, 0,
              -Y, -X, -(X+Y)*mu,  mu,  mu, -1,
              -Y,  X, -(X+Y)*mu,  mu, -mu, -1,
              Y, -X, -(X+Y)*mu, -mu,  mu, -1,
              Y,  X, -(X+Y)*mu, -mu, -mu, -1,
              Y,  X, -(X+Y)*mu,  mu,  mu,  1,
              Y, -X, -(X+Y)*mu,  mu, -mu,  1,
              -Y,  X, -(X+Y)*mu, -mu,  mu,  1,
              -Y, -X, -(X+Y)*mu, -mu, -mu,  1;
  }

  virtual void TearDown() {
  }

  void test_kinematics(Robot& robot, const ContactStatus& contact_status) const;
  void test_isFeasible(Robot& robot, const ContactStatus& contact_status) const;
  void test_setSlack(Robot& robot, const ContactStatus& contact_status) const;
  void test_evalConstraint(Robot& robot, const ContactStatus& contact_status) const;
  void test_evalDerivatives(Robot& robot, const ContactStatus& contact_status) const;
  void test_condenseSlackAndDual(Robot& robot, 
                                const ContactStatus& contact_status) const;
  void test_expandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const;

  double barrier, dt, mu, X, Y, fraction_to_boundary_rule;
  Eigen::MatrixXd cone;
};


void WrenchFrictionConeTest::test_kinematics(Robot& robot, 
                                             const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  EXPECT_TRUE(constr.useKinematics());
  EXPECT_TRUE(constr.kinematicsLevel() == KinematicsLevel::AccelerationLevel);
}


void WrenchFrictionConeTest::test_isFeasible(Robot& robot, 
                                             const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  EXPECT_EQ(constr.dimc(), 17*robot.maxNumSurfaceContacts());
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  bool feasible = true;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          const Eigen::VectorXd res = cone * s.f[i];
          if (res.maxCoeff() > 0) {
            feasible = false;
          }
        }
        break;
      default:
        break;
    }
  }
  EXPECT_EQ(constr.isFeasible(robot, contact_status, data, s), feasible);
}


void WrenchFrictionConeTest::test_setSlack(Robot& robot, const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.barrier()), data_ref(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  constr.allocateExtraData(data_ref);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateFrameKinematics(s.q);
  constr.setSlack(robot, contact_status, data, s);
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          data_ref.residual.segment(c_begin, 17).noalias() = cone * s.f[i];
          data_ref.slack.segment(c_begin, 17) = - data_ref.residual.segment(c_begin, 17);
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void WrenchFrictionConeTest::test_evalConstraint(Robot& robot, const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  ConstraintComponentData data(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  constr.evalConstraint(robot, contact_status, data, s);
  data_ref.residual.setZero();
  data_ref.cmpl.setZero();
  data_ref.log_barrier = 0;
  int c_begin = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          data_ref.residual.segment(c_begin, 17).noalias() 
              = cone * s.f[i] + data_ref.slack.segment(c_begin, 17);
          pdipm::computeComplementarySlackness(barrier, data_ref, c_begin, 17);
          data_ref.log_barrier += pdipm::logBarrier(barrier, data_ref.slack.segment(c_begin, 17));
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void WrenchFrictionConeTest::test_evalDerivatives(Robot& robot, const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, contact_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  constr.evalConstraint(robot, contact_status, data, s);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, contact_status, data, s, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          kkt_res_ref.lf().segment(dimf_stack, 6).noalias()
              += cone.transpose() * data_ref.dual.segment(c_begin, 17);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
}


void WrenchFrictionConeTest::test_condenseSlackAndDual(Robot& robot, 
                                                       const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, contact_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto kkt_mat = SplitKKTMatrix::Random(robot, contact_status);
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  constr.evalConstraint(robot, contact_status, data, s);
  constr.evalDerivatives(robot, contact_status, data, s, kkt_res);
  auto data_ref = data;
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(contact_status, data, kkt_mat, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          data_ref.r[0].array() = data_ref.dual.segment(c_begin, 17).array() 
                                    / data_ref.slack.segment(c_begin, 17).array();
          kkt_mat_ref.Qff().block(dimf_stack, dimf_stack, 6, 6).noalias()
              += cone.transpose() * data_ref.r[0].asDiagonal() * cone;
          pdipm::computeCondensingCoeffcient(data_ref, c_begin, 17);
          kkt_res_ref.lf().segment(dimf_stack, 6).noalias()
              += cone.transpose() * data_ref.cond.segment(c_begin, 17);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(kkt_res.isApprox(kkt_res_ref));
  EXPECT_TRUE(kkt_mat.isApprox(kkt_mat_ref));
}


void WrenchFrictionConeTest::test_expandSlackAndDual(Robot& robot, const ContactStatus& contact_status) const {
  WrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.barrier());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, contact_status);
  constr.setSlack(robot, contact_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  data.dslack.setRandom();
  data.ddual.setRandom();
  auto kkt_mat = SplitKKTMatrix::Random(robot, contact_status);
  auto kkt_res = SplitKKTResidual::Random(robot, contact_status);
  constr.evalConstraint(robot, contact_status, data, s);
  constr.evalDerivatives(robot, contact_status, data, s, kkt_res);
  constr.condenseSlackAndDual(contact_status, data, kkt_mat, kkt_res);
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot, contact_status);
  constr.expandSlackAndDual(contact_status, data, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<contact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (contact_status.isContactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (contact_status.isContactActive(i)) {
          data_ref.dslack.segment(c_begin, 17).noalias()
              = - cone * d.df().segment(dimf_stack, 6) 
                - data_ref.residual.segment(c_begin, 17);
          pdipm::computeDualDirection(data_ref, c_begin, 17);
          dimf_stack += 6;
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


TEST_F(WrenchFrictionConeTest, testWithHumanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto contact_status = robot.createContactStatus();
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
  contact_status.activateContact(1);
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
  contact_status.activateContacts({0, 1});
  test_kinematics(robot, contact_status);
  test_isFeasible(robot, contact_status);
  test_setSlack(robot, contact_status);
  test_evalConstraint(robot, contact_status);
  test_evalDerivatives(robot, contact_status);
  test_condenseSlackAndDual(robot, contact_status);
  test_expandSlackAndDual(robot, contact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}