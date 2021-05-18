#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/unconstrained_dynamics.hpp"

#include "robot_factory.hpp"


namespace idocp {

class UnconstrainedDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    robot = testhelper::CreateFixedBaseRobot();
    dt = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  Robot robot;
  double dt;
};


TEST_F(UnconstrainedDynamicsTest, linearizeUnconstrainedDynamics) {
  const SplitSolution s = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom(); 
  kkt_residual.la.setRandom();   
  kkt_residual.lu().setRandom(); 
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  UnconstrainedDynamics ud(robot);
  ud.linearizeUnconstrainedDynamics(robot, dt, s, kkt_residual);
  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  kkt_residual_ref.lq() += dt * dIDdq.transpose() * s.beta;
  kkt_residual_ref.lv() += dt * dIDdv.transpose() * s.beta;
  kkt_residual_ref.la   += dt * dIDda.transpose() * s.beta;
  kkt_residual_ref.lu() -= dt * s.beta;
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  EXPECT_DOUBLE_EQ(ud.l1NormUnconstrainedDynamicsResidual(dt), 
                   dt*ID.lpNorm<1>());
  EXPECT_DOUBLE_EQ(ud.squaredNormUnconstrainedDynamicsResidual(dt), 
                   dt*dt*ID.squaredNorm());
}


TEST_F(UnconstrainedDynamicsTest, condenseUnconstrainedDynamics) {
  const SplitSolution s = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom(); 
  kkt_residual.la.setRandom();   
  kkt_residual.lu().setRandom(); 
  UnconstrainedDynamics ud(robot);
  ud.linearizeUnconstrainedDynamics(robot, dt, s, kkt_residual);
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq().setRandom();
  kkt_matrix.Qvv().diagonal().setRandom();
  kkt_matrix.Qaa().diagonal().setRandom();
  kkt_matrix.Quu().diagonal().setRandom();

  SplitUnKKTResidual unkkt_residual(robot); 
  unkkt_residual.KKT_residual.setRandom();
  SplitUnKKTResidual unkkt_residual_ref = unkkt_residual;
  SplitUnKKTMatrix unkkt_matrix(robot);
  unkkt_matrix.Q.setRandom();
  SplitUnKKTMatrix unkkt_matrix_ref = unkkt_matrix;
  ud.condenseUnconstrainedDynamics(kkt_matrix, kkt_residual, 
                                   unkkt_matrix, unkkt_residual);

  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  const Eigen::VectorXd lu = kkt_residual.lu();
  const Eigen::VectorXd lu_condensed = kkt_residual.lu() + kkt_matrix.Quu() * ID;
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  unkkt_residual_ref.lq() = kkt_residual.lq() + dIDdq.transpose() * lu_condensed;
  unkkt_residual_ref.lv() = kkt_residual.lv() + dIDdv.transpose() * lu_condensed;
  unkkt_residual_ref.la() = kkt_residual.la   + dIDda.transpose() * lu_condensed;
  unkkt_residual_ref.Fx() = kkt_residual.Fx();
  unkkt_matrix_ref.Qqq() = kkt_matrix.Qqq();
  unkkt_matrix_ref.Qvv() = kkt_matrix.Qvv();
  unkkt_matrix_ref.Qaa() = kkt_matrix.Qaa();
  const Eigen::MatrixXd Quu_dIDdq = kkt_matrix.Quu() * dIDdq;
  const Eigen::MatrixXd Quu_dIDdv = kkt_matrix.Quu() * dIDdv;
  const Eigen::MatrixXd Quu_dIDda = kkt_matrix.Quu() * dIDda;
  unkkt_matrix_ref.Qqq() += dIDdq.transpose() * Quu_dIDdq;
  unkkt_matrix_ref.Qqv() =  dIDdq.transpose() * Quu_dIDdv;
  unkkt_matrix_ref.Qvv() += dIDdv.transpose() * Quu_dIDdv;
  unkkt_matrix_ref.Qaq() =  dIDda.transpose() * Quu_dIDdq;
  unkkt_matrix_ref.Qav() =  dIDda.transpose() * Quu_dIDdv;
  unkkt_matrix_ref.Qaa() += dIDda.transpose() * Quu_dIDda;
  EXPECT_TRUE(unkkt_residual_ref.isApprox(unkkt_residual));
  EXPECT_TRUE(unkkt_matrix_ref.isApprox(unkkt_matrix));

  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  ud.computeCondensedDirection(dt, kkt_matrix, kkt_residual, d);
  d_ref.du = ID + dIDdq * d_ref.dq() + dIDdv * d_ref.dv() + dIDda * d_ref.da();
  d_ref.dbeta() = (kkt_residual.lu() + kkt_matrix.Quu() * d_ref.du) / dt;
  EXPECT_TRUE(d_ref.isApprox(d));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}