#include <gtest/gtest.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/dynamics/terminal_state_equation.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class TerminalStateEquationTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  double dt;
};



TEST_F(TerminalStateEquationTest, fixedbase) {
  auto robot = testhelper::CreateRobotManipulator(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, q_prev, s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isZero());
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  state_equation.correctLinearizedStateEquation(kkt_matrix);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(TerminalStateEquationTest, floatingBase) {
  auto robot = testhelper::CreateQuadrupedalRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  TerminalStateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, q_prev, s, kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractConfiguration_dq0(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isZero());
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  state_equation.correctLinearizedStateEquation(kkt_matrix);
  const Eigen::MatrixXd dsubtract_dq_prev_inv = dsubtract_dq_prev.inverse();
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