#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/ocp/state_equation.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class StateEquationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dt;
};


TEST_F(StateEquationTest, fixedbase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  StateEquation state_equation(robot);
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  state_equation.linearizeStateEquation(robot, dt, q_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  kkt_residual_ref.Fq() = s.q + dt * s.v - s_next.q;
  kkt_residual_ref.Fv() = s.v + dt * s.a - s_next.v;
  kkt_residual_ref.lq() = s_next.lmd - s.lmd;
  kkt_residual_ref.lv() = dt * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual_ref.la   = dt * s_next.gmm;
  kkt_residual_ref.h  = s_next.lmd.dot(s.v) + s_next.gmm.dot(s.a);
  kkt_matrix_ref.hv() = s_next.lmd;
  kkt_matrix_ref.ha   = s_next.gmm;
  kkt_matrix_ref.fq() = s.v;
  kkt_matrix_ref.fv() = s.a;
  kkt_matrix_ref.Fqq() = Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fqv() = dt * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  state_equation.correctLinearizedStateEquation(robot, dt, s, s_next, 
                                                kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(StateEquationTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  StateEquation state_equation(robot);
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  state_equation.linearizeStateEquation(robot, dt, q_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractConfiguration_dqf(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractConfiguration_dq0(q_prev, s.q, dsubtract_dq_prev);
  kkt_residual_ref.Fq() = qdiff + dt * s.v;
  kkt_residual_ref.Fv() = s.v + dt * s.a - s_next.v;
  kkt_residual_ref.lq() = dsubtract_dq.transpose() * s_next.lmd
                          + dsubtract_dq_prev.transpose() * s.lmd;
  kkt_residual_ref.lv() = dt * s_next.lmd + s_next.gmm - s.gmm;
  kkt_residual_ref.la   = dt * s_next.gmm;
  kkt_residual_ref.h  = s_next.lmd.dot(s.v) + s_next.gmm.dot(s.a);
  kkt_matrix_ref.hv() = s_next.lmd;
  kkt_matrix_ref.ha   = s_next.gmm;
  kkt_matrix_ref.fq() = s.v;
  kkt_matrix_ref.fv() = s.a;
  kkt_matrix_ref.Fqq() = dsubtract_dq;
  kkt_matrix_ref.Fqv() = dt * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  kkt_matrix_ref.Fqq_prev = dsubtract_dq_prev;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  state_equation.correctLinearizedStateEquation(robot, dt, s, s_next, 
                                                kkt_matrix, kkt_residual);
  const Eigen::MatrixXd dsubtract_dq_prev_inv = dsubtract_dq_prev.topLeftCorner(6, 6).inverse();
  robot.dSubtractConfiguration_dq0(s.q, s_next.q, dsubtract_dq_prev);
  const Eigen::MatrixXd dsubtract_dq_inv = dsubtract_dq_prev.topLeftCorner(6, 6).inverse();
  Eigen::MatrixXd Fqq_ref = dsubtract_dq;
  Fqq_ref.topLeftCorner(6, 6) 
      = - dsubtract_dq_inv.topLeftCorner(6, 6) * dsubtract_dq.topLeftCorner(6, 6);
  Eigen::VectorXd Fq_ref = qdiff+dt*s.v;
  Eigen::MatrixXd Fqv_ref = dt * Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv());
  Fqv_ref.topLeftCorner(6, 6) = - dt * dsubtract_dq_inv.topLeftCorner(6, 6);
  Fq_ref.head(6) = - dsubtract_dq_inv.topLeftCorner(6, 6) * (qdiff+dt*s.v).head(6);
  kkt_matrix_ref.Fqq() = Fqq_ref;
  kkt_matrix_ref.Fqv() = Fqv_ref;
  kkt_matrix_ref.Fqq_prev = dsubtract_dq_prev;
  kkt_residual_ref.Fq() = Fq_ref;
  const Eigen::VectorXd fq_tmp = kkt_matrix_ref.fq().head(6);
  kkt_matrix_ref.fq().head(6) = - dsubtract_dq_inv.topLeftCorner(6, 6) * fq_tmp;
  EXPECT_TRUE(kkt_residual.isApprox(kkt_residual_ref));
  EXPECT_TRUE(kkt_matrix.isApprox(kkt_matrix_ref));
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  d_ref.dlmd().head(6) 
      = - dsubtract_dq_prev_inv.topLeftCorner(6, 6).transpose() * d_ref.dlmd().head(6);
  EXPECT_TRUE(d.isApprox(d_ref));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}