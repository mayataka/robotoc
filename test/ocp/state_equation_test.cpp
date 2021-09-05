#include <gtest/gtest.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/state_equation.hpp"

#include "robot_factory.hpp"


namespace idocp {

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
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dt*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  kkt_matrix.setZero();
  kkt_residual.setZero();
  state_equation.linearizeStateEquationAlongLieGroup(robot, dt, q_prev, s, s_next, 
                                                     kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dt*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  state_equation.correctCostateDirection(d);
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(StateEquationTest, floatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  const Eigen::VectorXd q_prev = robot.generateFeasibleConfiguration();
  const auto s = SplitSolution::Random(robot);
  const auto s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  StateEquation state_equation(robot);
  state_equation.linearizeStateEquation(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractConfiguration_dqf(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractConfiguration_dq0(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  kkt_matrix.setZero();
  kkt_residual.setZero();
  state_equation.linearizeStateEquationAlongLieGroup(robot, dt, q_prev, s, s_next, 
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
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Fqq_ref));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(Fqv_ref));
  EXPECT_TRUE(kkt_residual.Fq().isApprox((Fq_ref)));
  auto d = SplitDirection::Random(robot);
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