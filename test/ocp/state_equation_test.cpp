#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
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


TEST_F(StateEquationTest, forwardEulerFixedbase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dt*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  stateequation::condenseForwardEuler(robot, dt, s, s_next.q, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((s.q+dt*s.v-s_next.q)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(StateEquationTest, forwardEulerTerminalFixedbase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeForwardEulerTerminal(robot, q_prev, s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isZero());
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  stateequation::condenseForwardEulerTerminal(robot, kkt_matrix);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(StateEquationTest, forwardEulerFloatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(s.q, s_next.q, qdiff);
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dq);
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((s.v+dt*s.a-s_next.v)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq.transpose()*s_next.lmd+dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s_next.lmd+s_next.gmm-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s_next.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dq));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(dt*Eigen::MatrixXd::Identity(robot.dimv(), robot.dimv())));
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  std::cout << "kkt_matrix.Fqq()" << std::endl;
  std::cout << kkt_matrix.Fqq() << std::endl;
  std::cout << "kkt_matrix.Fqq_prev" << std::endl;
  std::cout << kkt_matrix.Fqq_prev << std::endl;
  stateequation::condenseForwardEuler(robot, dt, s, s_next.q, kkt_matrix, kkt_residual);
  Eigen::MatrixXd dsubtract_dq_prev_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dq_prev, dsubtract_dq_prev_inv);
  robot.dSubtractdConfigurationMinus(s.q, s_next.q, dsubtract_dq_prev);
  Eigen::MatrixXd dsubtract_dq_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dq_prev, dsubtract_dq_inv);
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
  EXPECT_TRUE(kkt_matrix.Fqq_inv.isApprox(dsubtract_dq_inv.topLeftCorner(6, 6)));
  EXPECT_TRUE(kkt_matrix.Fqq_prev_inv.isApprox(dsubtract_dq_prev_inv.topLeftCorner(6, 6)));
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  dlmd_ref.head(6) = - dsubtract_dq_prev_inv.topLeftCorner(6, 6).transpose() * dlmd_ref.head(6);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(StateEquationTest, forwardEulerTerminalFloatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeForwardEulerTerminal(robot, q_prev, s, kkt_matrix, kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  Eigen::MatrixXd dsubtract_dq = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dq_prev = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dq_prev);
  EXPECT_TRUE(kkt_residual.Fq().isZero());
  EXPECT_TRUE(kkt_residual.Fv().isZero());
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dq_prev.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isZero());
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_TRUE(kkt_matrix.Fqv().isZero());
  EXPECT_TRUE(kkt_matrix.Fqq_prev.isApprox(dsubtract_dq_prev));
  stateequation::condenseForwardEulerTerminal(robot, kkt_matrix);
  Eigen::MatrixXd dsubtract_dq_prev_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dq_prev, dsubtract_dq_prev_inv);
  EXPECT_TRUE(kkt_matrix.Fqq_prev_inv.isApprox(dsubtract_dq_prev_inv.topLeftCorner(6, 6)));
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionForwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  dlmd_ref.head(6) = - dsubtract_dq_prev_inv.topLeftCorner(6, 6).transpose() * dlmd_ref.head(6);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(StateEquationTest, backwardEulerFixedBase) {
  auto robot = testhelper::CreateFixedBaseRobot(dt);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dt*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((s_next.lmd-s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::linearizeBackwardEulerTerminal(robot, dt, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dt*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  stateequation::condenseBackwardEuler(robot, dt, q_prev, s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((q_prev-s.q+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dt*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((-1*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isZero());
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionBackwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}


TEST_F(StateEquationTest, backwardEulerFloatingBase) {
  auto robot = testhelper::CreateFloatingBaseRobot(dt);
  SplitSolution s = SplitSolution::Random(robot);
  SplitSolution s_next = SplitSolution::Random(robot);
  Eigen::VectorXd q_prev = Eigen::VectorXd::Random(robot.dimq());
  robot.normalizeConfiguration(q_prev);
  const Eigen::VectorXd v_prev = Eigen::VectorXd::Random(robot.dimv());
  SplitKKTResidual kkt_residual(robot);
  SplitKKTMatrix kkt_matrix(robot);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next,
                                        kkt_matrix,  kkt_residual);
  Eigen::VectorXd qdiff = Eigen::VectorXd::Zero(robot.dimv());
  robot.subtractConfiguration(q_prev, s.q, qdiff);
  Eigen::MatrixXd dsubtract_dqminus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Eigen::MatrixXd dsubtract_dqplus = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationMinus(q_prev, s.q, dsubtract_dqminus);
  robot.dSubtractdConfigurationPlus(s.q, s_next.q, dsubtract_dqplus);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dt*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqplus.transpose()*s_next.lmd+dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s.lmd-s.gmm+s_next.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  kkt_residual.setZero();
  kkt_matrix.setZero();
  stateequation::linearizeBackwardEulerTerminal(robot, dt, q_prev, v_prev, 
                                                s, kkt_matrix, kkt_residual);
  EXPECT_TRUE(kkt_residual.Fq().isApprox((qdiff+dt*s.v)));
  EXPECT_TRUE(kkt_residual.Fv().isApprox((v_prev-s.v+dt*s.a)));
  EXPECT_TRUE(kkt_residual.lq().isApprox((dsubtract_dqminus.transpose()*s.lmd)));
  EXPECT_TRUE(kkt_residual.lv().isApprox((dt*s.lmd-s.gmm)));
  EXPECT_TRUE(kkt_residual.la.isApprox((dt*s.gmm)));
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(dsubtract_dqminus));
  std::cout << "kkt_matrix.Fqq()" << std::endl;
  std::cout << kkt_matrix.Fqq() << std::endl;
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.lpNorm<1>(), 
                   stateequation::l1NormStateEuqationResidual(kkt_residual));
  EXPECT_DOUBLE_EQ(kkt_residual.Fx.squaredNorm(), 
                   stateequation::squaredNormStateEuqationResidual(kkt_residual));
  stateequation::condenseBackwardEuler(robot, dt, q_prev, s, kkt_matrix, kkt_residual);
  robot.dSubtractdConfigurationPlus(q_prev, s.q, dsubtract_dqplus);
  Eigen::MatrixXd dsubtract_dq_inv = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  robot.dSubtractdConfigurationInverse(dsubtract_dqplus, dsubtract_dq_inv);
  Eigen::MatrixXd Fqq_ref = dsubtract_dqminus;
  Fqq_ref.topLeftCorner(6, 6) 
      = dsubtract_dq_inv.topLeftCorner(6, 6) * dsubtract_dqminus.topLeftCorner(6, 6);
  Eigen::VectorXd Fq_ref = qdiff+dt*s.v;
  Eigen::MatrixXd Fqv_ref = Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv());
  Fqv_ref.topLeftCorner(6, 6) = dt * dsubtract_dq_inv.topLeftCorner(6, 6);
  Fq_ref.head(6) = dsubtract_dq_inv.topLeftCorner(6, 6) * (qdiff+dt*s.v).head(6);
  EXPECT_TRUE(kkt_matrix.Fqq().isApprox(Fqq_ref));
  EXPECT_TRUE(kkt_matrix.Fqv().isApprox(Fqv_ref));
  EXPECT_TRUE(kkt_residual.Fq().isApprox((Fq_ref)));
  EXPECT_TRUE(kkt_matrix.Fqq_inv.isApprox(dsubtract_dq_inv.topLeftCorner(6, 6)));
  Eigen::VectorXd dlmd = Eigen::VectorXd::Random(robot.dimv());
  Eigen::VectorXd dlmd_ref = dlmd;
  stateequation::correctCostateDirectionBackwardEuler(robot, kkt_matrix, kkt_residual, dlmd);
  dlmd_ref.head(6) = dsubtract_dq_inv.topLeftCorner(6, 6).transpose() * dlmd_ref.head(6);
  EXPECT_TRUE(dlmd.isApprox(dlmd_ref));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}