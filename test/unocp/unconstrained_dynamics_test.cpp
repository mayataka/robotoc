#include <string>
#include <iostream>

#include <gtest/gtest.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/unocp/unconstrained_dynamics.hpp"


namespace idocp {

class UnconstrainedDynamicsTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf = "../urdf/iiwa14/iiwa14.urdf";
    robot = Robot(urdf);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  std::string urdf;
  Robot robot;
  double dtau;
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
  ud.linearizeUnconstrainedDynamics(robot, dtau, s, kkt_residual);
  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  kkt_residual_ref.lq() += dtau * dIDdq.transpose() * s.beta;
  kkt_residual_ref.lv() += dtau * dIDdv.transpose() * s.beta;
  kkt_residual_ref.la   += dtau * dIDda.transpose() * s.beta;
  kkt_residual_ref.lu() -= dtau * s.beta;
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  EXPECT_DOUBLE_EQ(ud.l1NormUnconstrainedDynamicsResidual(dtau), 
                   dtau*ID.lpNorm<1>());
  EXPECT_DOUBLE_EQ(ud.squaredNormUnconstrainedDynamicsResidual(dtau), 
                   dtau*dtau*ID.squaredNorm());
}


TEST_F(UnconstrainedDynamicsTest, condenseUnconstrainedDynamics) {
  const SplitSolution s = SplitSolution::Random(robot);
  SplitKKTResidual kkt_residual(robot);
  kkt_residual.lq().setRandom();
  kkt_residual.lv().setRandom(); 
  kkt_residual.la.setRandom();   
  kkt_residual.lu().setRandom(); 
  UnconstrainedDynamics ud(robot);
  ud.linearizeUnconstrainedDynamics(robot, dtau, s, kkt_residual);
  SplitKKTResidual kkt_residual_ref = kkt_residual;
  SplitKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq().setRandom();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq().setRandom();
  kkt_matrix.Qvv().setRandom();
  kkt_matrix.Qaa().setRandom();
  kkt_matrix.Quu().diagonal().setRandom();
  SplitKKTMatrix kkt_matrix_ref = kkt_matrix;
  Eigen::MatrixXd Qaq(Eigen::MatrixXd::Random(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd Qav(Eigen::MatrixXd::Random(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd Qaq_ref = Qaq;
  Eigen::MatrixXd Qav_ref = Qav;
  ud.condenseUnconstrainedDynamics(robot, dtau, kkt_matrix, kkt_residual, Qaq, Qav);

  Eigen::VectorXd ID(Eigen::VectorXd::Zero(robot.dimv()));
  robot.RNEA(s.q, s.v, s.a, ID);
  ID.noalias() -= s.u;
  const Eigen::VectorXd lu_origin = kkt_residual_ref.lu();
  const Eigen::VectorXd lu_condensed = kkt_residual_ref.lu() + kkt_matrix_ref.Quu() * ID;
  Eigen::MatrixXd dIDdq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  Eigen::MatrixXd dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv()));
  robot.RNEADerivatives(s.q, s.v, s.a, dIDdq, dIDdv, dIDda);
  kkt_residual_ref.lq().noalias() += dIDdq.transpose() * lu_condensed;
  kkt_residual_ref.lv().noalias() += dIDdv.transpose() * lu_condensed;
  kkt_residual_ref.la.noalias() += dIDda.transpose() * lu_condensed;

  const Eigen::MatrixXd Quu_dIDda = kkt_matrix_ref.Quu().diagonal().asDiagonal() * dIDda;
  kkt_matrix_ref.Qaa().noalias() += dIDda.transpose() * Quu_dIDda;
  const Eigen::MatrixXd Quu_dIDdq = kkt_matrix_ref.Quu().diagonal().asDiagonal() * dIDdq;
  const Eigen::MatrixXd Quu_dIDdv = kkt_matrix_ref.Quu().diagonal().asDiagonal() * dIDdv;
  Qaq_ref = dIDda.transpose() * Quu_dIDdq;
  Qav_ref = dIDda.transpose() * Quu_dIDdv;
  kkt_matrix_ref.Qqq().noalias() += dIDdq.transpose() * Quu_dIDdq;
  kkt_matrix_ref.Qqv().noalias() += dIDdq.transpose() * Quu_dIDdv;
  kkt_matrix_ref.Qvv().noalias() += dIDdv.transpose() * Quu_dIDdv;
  EXPECT_TRUE(kkt_residual_ref.isApprox(kkt_residual));
  EXPECT_TRUE(kkt_matrix_ref.isApprox(kkt_matrix));

  SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_ref = d;
  ud.computeCondensedDirection(dtau, kkt_matrix, kkt_residual, d);
  d_ref.du() = ID + dIDdq * d_ref.dq() + dIDdv * d_ref.dv() + dIDda * d_ref.da();
  d_ref.dbeta() = (kkt_residual_ref.lu() + kkt_matrix_ref.Quu() * d_ref.du()) / dtau;
  EXPECT_TRUE(d_ref.isApprox(d));
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}