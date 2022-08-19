#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"


#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/constraints/pdipm.hpp"
#include "robotoc/constraints/impulse_wrench_friction_cone.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class ImpulseWrenchFrictionConeTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    barrier_param = 1.0e-03;
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

  void test_isFeasible(Robot& robot, const ImpulseStatus& impulse_status) const;
  void test_setSlack(Robot& robot, const ImpulseStatus& impulse_status) const;
  void test_evalConstraint(Robot& robot, const ImpulseStatus& impulse_status) const;
  void test_evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status) const;
  void test_condenseSlackAndDual(Robot& robot, 
                                const ImpulseStatus& impulse_status) const;
  void test_expandSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const;

  double barrier_param, dt, mu, X, Y, fraction_to_boundary_rule;
  Eigen::MatrixXd cone;
};


void ImpulseWrenchFrictionConeTest::test_isFeasible(Robot& robot, 
                                                    const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  EXPECT_EQ(constr.dimc(), 17*robot.maxNumSurfaceContacts());
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateFrameKinematics(s.q);
  bool feasible = true;
  for (int i=0; i<impulse_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
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
  EXPECT_EQ(constr.isFeasible(robot, impulse_status, data, s), feasible);
}


void ImpulseWrenchFrictionConeTest::test_setSlack(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam()), data_ref(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  constr.allocateExtraData(data_ref);
  const int dimc = constr.dimc();
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateFrameKinematics(s.q);
  constr.setSlack(robot, impulse_status, data, s);
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
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


void ImpulseWrenchFrictionConeTest::test_evalConstraint(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  const int dimc = constr.dimc();
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q);
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto data_ref = data;
  constr.evalConstraint(robot, impulse_status, data, s);
  data_ref.residual.setZero();
  data_ref.cmpl.setZero();
  data_ref.log_barrier = 0;
  int c_begin = 0;
  for (int i=0; i<impulse_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
          data_ref.residual.segment(c_begin, 17).noalias() 
              = cone * s.f[i] + data_ref.slack.segment(c_begin, 17);
          pdipm::computeComplementarySlackness(barrier_param, data_ref, c_begin, 17);
          data_ref.log_barrier += pdipm::logBarrier(barrier_param, data_ref.slack.segment(c_begin, 17));
        }
        c_begin += 17;
        break;
      default:
        break;
    }
  }
  EXPECT_TRUE(data.isApprox(data_ref));
}


void ImpulseWrenchFrictionConeTest::test_evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, impulse_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  constr.evalConstraint(robot, impulse_status, data, s);
  auto data_ref = data;
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot, impulse_status);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, impulse_status, data, s, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
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


void ImpulseWrenchFrictionConeTest::test_condenseSlackAndDual(Robot& robot, 
                                                              const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, impulse_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot, impulse_status);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot, impulse_status);
  constr.evalConstraint(robot, impulse_status, data, s);
  constr.evalDerivatives(robot, impulse_status, data, s, kkt_res);
  auto data_ref = data;
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(impulse_status, data, kkt_mat, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
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


void ImpulseWrenchFrictionConeTest::test_expandSlackAndDual(Robot& robot, const ImpulseStatus& impulse_status) const {
  ImpulseWrenchFrictionCone constr(robot, mu, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = ImpulseSplitSolution::Random(robot, impulse_status);
  constr.setSlack(robot, impulse_status, data, s);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  data.dslack.setRandom();
  data.ddual.setRandom();
  auto kkt_mat = ImpulseSplitKKTMatrix::Random(robot, impulse_status);
  auto kkt_res = ImpulseSplitKKTResidual::Random(robot, impulse_status);
  constr.evalConstraint(robot, impulse_status, data, s);
  constr.evalDerivatives(robot, impulse_status, data, s, kkt_res);
  constr.condenseSlackAndDual(impulse_status, data, kkt_mat, kkt_res);
  auto data_ref = data;
  const auto d = ImpulseSplitDirection::Random(robot, impulse_status);
  constr.expandSlackAndDual(impulse_status, data, d);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<impulse_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impulse_status.isImpulseActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impulse_status.isImpulseActive(i)) {
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


TEST_F(ImpulseWrenchFrictionConeTest, testWithHumanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto impulse_status = robot.createImpulseStatus();
  test_isFeasible(robot, impulse_status);
  test_setSlack(robot, impulse_status);
  test_evalConstraint(robot, impulse_status);
  test_evalDerivatives(robot, impulse_status);
  test_condenseSlackAndDual(robot, impulse_status);
  test_expandSlackAndDual(robot, impulse_status);
  impulse_status.activateImpulse(1);
  test_isFeasible(robot, impulse_status);
  test_setSlack(robot, impulse_status);
  test_evalConstraint(robot, impulse_status);
  test_evalDerivatives(robot, impulse_status);
  test_condenseSlackAndDual(robot, impulse_status);
  test_expandSlackAndDual(robot, impulse_status);
  impulse_status.activateImpulses({0, 1});
  test_isFeasible(robot, impulse_status);
  test_setSlack(robot, impulse_status);
  test_evalConstraint(robot, impulse_status);
  test_evalDerivatives(robot, impulse_status);
  test_condenseSlackAndDual(robot, impulse_status);
  test_expandSlackAndDual(robot, impulse_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}