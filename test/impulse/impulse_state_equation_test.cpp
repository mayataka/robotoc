#include <gtest/gtest.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
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


TEST_F(ImpulseStateEquationTest, fixedBase) {
  const double dt = 0.001;
  const auto robot = testhelper::CreateFixedBaseRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = ImpulseSplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isIdentity());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  state_equation.linearizeStateEquationAlongLieGroup(robot, q_prev, s, s_next, 
                                                     kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isIdentity());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  auto d = ImpulseSplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ImpulseStateEquationTest, floatingBase) {
  const double dt = 0.001;
  const auto robot = testhelper::CreateFloatingBaseRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const ImpulseSplitSolution s = ImpulseSplitSolution::Random(robot);
  const SplitSolution s_next = SplitSolution::Random(robot);
  ImpulseSplitKKTResidual kkt_residual(robot);
  ImpulseSplitKKTMatrix kkt_matrix(robot);
  ImpulseStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, q_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractConfiguration_dqf(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractConfiguration_dq0(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+s.dv-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.ldv.isApprox((s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  kkt_matrix.setZero();
  kkt_residual.setZero();
  state_equation.linearizeStateEquationAlongLieGroup(robot, q_prev, s, s_next, 
                                                     kkt_matrix, kkt_residual);
  const Eigen::MatrixXd dsubtract_dq_prev_inv = dsubtract_dq_prev.topLeftCorner(6, 6).inverse();
  robot.dSubtractConfiguration_dq0(s.q, s_next.q, dsubtract_dq_prev);
  Eigen::MatrixXd dsubtract_dq_inv = dsubtract_dq_prev.topLeftCorner(6, 6).inverse();
  Eigen::MatrixXd Fqq_ref = dsubtract_dq;
  Fqq_ref.topLeftCorner(6, 6) 
      = - dsubtract_dq_inv.topLeftCorner(6, 6) * dsubtract_dq.topLeftCorner(6, 6);
  Eigen::VectorXd Fq_ref = qdiff;
  Eigen::MatrixXd Fqv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Fq_ref.head(6) = - dsubtract_dq_inv.topLeftCorner(6, 6) * qdiff.head(6);
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Fqq_ref));
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_residual.Fq().isApprox((Fq_ref)));
  auto d = ImpulseSplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  d_ref.dlmd().head(6) 
      = - dsubtract_dq_prev_inv.topLeftCorner(6, 6).transpose() * d_ref.dlmd().head(6);
  EXPECT_TRUE(d.isApprox(d_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}