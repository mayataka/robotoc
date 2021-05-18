#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_state_equation.hpp"

#include "robot_factory.hpp"


namespace idocp {

class ImpulseStateEquationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
  }

  virtual void TearDown() {
  }
};


TEST_F(ImpulseStateEquationTest, forwardEulerFixedBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  const auto s = ImpulseSplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  stateequation::condenseImpulseForwardEuler(robot, s, s_next.q, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  // stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  stateequation::correctCostateDirectionImpulseForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(ImpulseStateEquationTest, forwardEulerFloatingBase) {
  const double dt = 0.001;
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeImpulseForwardEuler(robot, q_prev, s, s_next, 
                                              kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  stateequation::condenseImpulseForwardEuler(robot, s, s_next.q, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dsubtract_dq_prev_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dq_prev, dsubtract_dq_prev_inv);
  robot.dSubtractdConfigurationMinus(s.q, s_next.q, dsubtract_dq_prev);
  Eigen::MatrixXd dsubtract_dq_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dq_prev, dsubtract_dq_inv);
  Eigen::MatrixXd Fqq_ref = dsubtract_dq;
  Fqq_ref.topLeftCorner(6, 6) 
      = - dsubtract_dq_inv.topLeftCorner(6, 6) * dsubtract_dq.topLeftCorner(6, 6);
  Eigen::VectorXd Fq_ref = qdiff;
  Eigen::MatrixXd Fqv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Fq_ref.head(6) = - dsubtract_dq_inv.topLeftCorner(6, 6) * qdiff.head(6);
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Fqq_ref));
  EXPECT_TRUE(kkt_residual.Fq().isApprox((Fq_ref)));
  EXPECT_TRUE(kkt_matrix.Fqq_inv.isApprox(dsubtract_dq_inv.topLeftCorner(robot.dim_passive(), robot.dim_passive())));
  EXPECT_TRUE(kkt_matrix.Fqq_prev_inv.isApprox(dsubtract_dq_prev_inv.topLeftCorner(robot.dim_passive(), robot.dim_passive())));
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  // stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  stateequation::correctCostateDirectionImpulseForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  dlmd_ref.head(6) = - dsubtract_dq_prev_inv.topLeftCorner(6, 6).transpose() * dlmd_ref.head(6);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}