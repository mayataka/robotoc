#include <gtest/gtest.h>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/dynamics/unconstr_dynamics.hpp"

#include "robot_factory.hpp"


namespace robotoc {

class UnconstrDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateRobotManipulator();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
};


TEST_F(UnconstrDynamicsTest, linearizeUnconstrDynamics) {
  const auto s = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx.setRandom();
  kkt_residual.la.setRandom();   
  kkt_residual.lu.setRandom(); 
  auto kkt_residual_ref = kkt_residual;
  UnconstrDynamics ud(robot);
  ud.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  kkt_residual_ref.lq() += dt * dIDdq.transpose() * s.beta;
  kkt_residual_ref.lv() += dt * dIDdv.transpose() * s.beta;
  kkt_residual_ref.la   += dt * dIDda.transpose() * s.beta;
  kkt_residual_ref.lu -= dt * s.beta;
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  EXPECT_DOUBLE_EQ(ud.primalFeasibility(), ID.lpNorm<1>());
  EXPECT_DOUBLE_EQ(ud.KKTError(), ID.squaredNorm());
}


TEST_F(UnconstrDynamicsTest, condenseUnconstrDynamics) {
  const SplitSolution s = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lx.setRandom();
  kkt_residual.la.setRandom();   
  kkt_residual.lu.setRandom(); 
  UnconstrDynamics ud(robot);
  ud.linearizeUnconstrDynamics(robot, dt, s, kkt_residual);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq().setRandom();
  kkt_matrix.Qvv().diagonal().setRandom();
  kkt_matrix.Qaa.diagonal().setRandom();
  kkt_matrix.Quu.diagonal().setRandom();
  auto kkt_residual_ref = kkt_residual;
  auto kkt_matrix_ref = kkt_matrix;
  ud.condenseUnconstrDynamics(kkt_matrix, kkt_residual);

  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  const Eigen::VectorXd lu = kkt_residual.lu;
  const Eigen::VectorXd lu_condensed = kkt_residual.lu + kkt_matrix.Quu * ID;
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  kkt_residual_ref.lq() += dIDdq.transpose() * lu_condensed;
  kkt_residual_ref.lv() += dIDdv.transpose() * lu_condensed;
  kkt_residual_ref.la   += dIDda.transpose() * lu_condensed;
  const Eigen::MatrixXd Quu_dIDdq = kkt_matrix.Quu * dIDdq;
  const Eigen::MatrixXd Quu_dIDdv = kkt_matrix.Quu * dIDdv;
  const Eigen::MatrixXd Quu_dIDda = kkt_matrix.Quu * dIDda;
  kkt_matrix_ref.Qqq() += dIDdq.transpose() * Quu_dIDdq;
  kkt_matrix_ref.Qqv() += dIDdq.transpose() * Quu_dIDdv;
  kkt_matrix_ref.Qvq()  = kkt_matrix_ref.Qqv().transpose();
  kkt_matrix_ref.Qvv() += dIDdv.transpose() * Quu_dIDdv;
  kkt_matrix_ref.Qqu().transpose() =  dIDda.transpose() * Quu_dIDdq;
  kkt_matrix_ref.Qvu().transpose() =  dIDda.transpose() * Quu_dIDdv;
  kkt_matrix_ref.Qaa   += dIDda.transpose() * Quu_dIDda;
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));

  auto d = SplitDirection::Random(robot);
  auto d_ref = d;
  ud.expandPrimal(d);
  d_ref.du = ID + dIDdq * d_ref.dq() + dIDdv * d_ref.dv() + dIDda * d_ref.da();
  EXPECT_TRUE(d_ref.isApprox(d));
  ud.expandDual(dt, kkt_matrix, kkt_residual, d);
  d_ref.dbeta() = (kkt_residual.lu + kkt_matrix.Quu * d_ref.du) / dt;
  EXPECT_TRUE(d_ref.isApprox(d));
}

} // namespace robotoc


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}