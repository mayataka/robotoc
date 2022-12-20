#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"


#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/constraints/pdipm.hpp"
#include "robotoc/constraints/impact_wrench_cone.hpp"

#include "robot_factory.hpp"

namespace robotoc {

class ImpactWrenchConeTest : public ::testing::Test {
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
    grid_info = GridInfo::Random();
    grid_info.type = GridType::Impact;
    grid_info.dt = 0.0;
  }

  virtual void TearDown() {
  }

  void test_isFeasible(Robot& robot, const ImpactStatus& impact_status) const;
  void test_setSlack(Robot& robot, const ImpactStatus& impact_status) const;
  void test_evalConstraint(Robot& robot, const ImpactStatus& impact_status) const;
  void test_evalDerivatives(Robot& robot, const ImpactStatus& impact_status) const;
  void test_condenseSlackAndDual(Robot& robot, 
                                const ImpactStatus& impact_status) const;
  void test_expandSlackAndDual(Robot& robot, const ImpactStatus& impact_status) const;

  double barrier_param, dt, mu, X, Y, fraction_to_boundary_rule;
  Eigen::MatrixXd cone;
  GridInfo grid_info;
};


void ImpactWrenchConeTest::test_isFeasible(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  EXPECT_EQ(constr.dimc(), 17*robot.maxNumSurfaceContacts());
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateFrameKinematics(s.q);
  bool feasible = true;
  for (int i=0; i<impact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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
  EXPECT_EQ(constr.isFeasible(robot, impact_status, grid_info, s, data), feasible);
}


void ImpactWrenchConeTest::test_setSlack(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam()), data_ref(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  constr.allocateExtraData(data_ref);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateFrameKinematics(s.q);
  constr.setSlack(robot, impact_status, grid_info, s, data);
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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


void ImpactWrenchConeTest::test_evalConstraint(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, impact_status);
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
  constr.evalConstraint(robot, impact_status, grid_info, s, data);
  data_ref.residual.setZero();
  data_ref.cmpl.setZero();
  data_ref.log_barrier = 0;
  int c_begin = 0;
  for (int i=0; i<impact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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


void ImpactWrenchConeTest::test_evalDerivatives(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, impact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  constr.evalConstraint(robot, impact_status, grid_info, s, data);
  auto data_ref = data;
  auto kkt_res = SplitKKTResidual::Random(robot, impact_status);
  auto kkt_res_ref = kkt_res;
  constr.evalDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impact_status.isImpactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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


void ImpactWrenchConeTest::test_condenseSlackAndDual(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, impact_status);
  robot.updateKinematics(s.q);
  constr.setSlack(robot, impact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  // auto kkt_mat = SplitKKTMatrix::Random(robot, impact_status);
  auto kkt_mat = SplitKKTMatrix(robot);
  kkt_mat.setContactDimension(impact_status.dimf());
  kkt_mat.setRandom();
  auto kkt_res = SplitKKTResidual::Random(robot, impact_status);
  constr.evalConstraint(robot, impact_status, grid_info, s, data);
  constr.evalDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  auto data_ref = data;
  auto kkt_mat_ref = kkt_mat;
  auto kkt_res_ref = kkt_res;
  constr.condenseSlackAndDual(impact_status, grid_info, data, kkt_mat, kkt_res);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<robot.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impact_status.isImpactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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


void ImpactWrenchConeTest::test_expandSlackAndDual(Robot& robot, const ImpactStatus& impact_status) const {
  ImpactWrenchCone constr(robot, X, Y); 
  ConstraintComponentData data(constr.dimc(), constr.getBarrierParam());
  constr.allocateExtraData(data);
  const int dimc = constr.dimc();
  const auto s = SplitSolution::Random(robot, impact_status);
  constr.setSlack(robot, impact_status, grid_info, s, data);
  data.slack.setRandom();
  data.dual.setRandom();
  data.slack = data.slack.array().abs();
  data.dual = data.dual.array().abs();
  data.residual.setRandom();
  data.cmpl.setRandom();
  data.dslack.setRandom();
  data.ddual.setRandom();
  // auto kkt_mat = SplitKKTMatrix::Random(robot, impact_status);
  auto kkt_mat = SplitKKTMatrix(robot);
  kkt_mat.setContactDimension(impact_status.dimf());
  kkt_mat.setRandom();
  auto kkt_res = SplitKKTResidual::Random(robot, impact_status);
  constr.evalConstraint(robot, impact_status, grid_info, s, data);
  constr.evalDerivatives(robot, impact_status, grid_info, s, data, kkt_res);
  constr.condenseSlackAndDual(impact_status, grid_info, data, kkt_mat, kkt_res);
  auto data_ref = data;
  const auto d = SplitDirection::Random(robot, impact_status);
  constr.expandSlackAndDual(impact_status, grid_info, d, data);
  data_ref.dslack.fill(1.0);
  data_ref.ddual.fill(1.0);
  int dimf_stack = 0;
  int c_begin = 0;
  for (int i=0; i<impact_status.maxNumContacts(); ++i) {
    switch (robot.contactType(i)) {
      case ContactType::PointContact:
        if (impact_status.isImpactActive(i)) {
          dimf_stack += 3;
        }
        break;
      case ContactType::SurfaceContact:
        if (impact_status.isImpactActive(i)) {
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


TEST_F(ImpactWrenchConeTest, testWithHumanoidRobot) {
  auto robot = testhelper::CreateHumanoidRobot(dt);
  auto impact_status = robot.createImpactStatus();
  test_isFeasible(robot, impact_status);
  test_setSlack(robot, impact_status);
  test_evalConstraint(robot, impact_status);
  test_evalDerivatives(robot, impact_status);
  test_condenseSlackAndDual(robot, impact_status);
  test_expandSlackAndDual(robot, impact_status);
  impact_status.activateImpact(1);
  test_isFeasible(robot, impact_status);
  test_setSlack(robot, impact_status);
  test_evalConstraint(robot, impact_status);
  test_evalDerivatives(robot, impact_status);
  test_condenseSlackAndDual(robot, impact_status);
  test_expandSlackAndDual(robot, impact_status);
  impact_status.activateImpacts({0, 1});
  test_isFeasible(robot, impact_status);
  test_setSlack(robot, impact_status);
  test_evalConstraint(robot, impact_status);
  test_evalDerivatives(robot, impact_status);
  test_condenseSlackAndDual(robot, impact_status);
  test_expandSlackAndDual(robot, impact_status);
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}