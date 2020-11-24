#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"


namespace idocp {

class ImpulseRiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
  }

  virtual void TearDown() {
  }

  static ImpulseKKTMatrix createKKTMatrix(const Robot& robot);
  static ImpulseKKTResidual createKKTResidual(const Robot& robot);
  static RiccatiFactorization createRiccatiFactorization(const Robot& robot);

  static void testBackwardRecursion(const Robot& robot);
  static void testForwardRecursion(const Robot& robot);
  static void testFactorizeStateConstraintFactorization(const Robot& robot);

  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


ImpulseKKTMatrix ImpulseRiccatiFactorizerTest::createKKTMatrix(const Robot& robot) {
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  ImpulseKKTMatrix kkt_matrix(robot);
  kkt_matrix.Qqq() = seed * seed.transpose();
  kkt_matrix.Qqv().setRandom();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  kkt_matrix.Qvv() = seed * seed.transpose();
  if (robot.has_floating_base()) {
    kkt_matrix.Fqq().setIdentity();
    kkt_matrix.Fqq().topLeftCorner(robot.dim_passive(), robot.dim_passive()).setRandom();
  }
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  return kkt_matrix;
}


ImpulseKKTResidual ImpulseRiccatiFactorizerTest::createKKTResidual(const Robot& robot) {
  ImpulseKKTResidual kkt_residual(robot);
  kkt_residual.lx().setRandom();
  kkt_residual.Fx().setRandom();
  return kkt_residual;
}


RiccatiFactorization ImpulseRiccatiFactorizerTest::createRiccatiFactorization(const Robot& robot) {
  RiccatiFactorization riccati(robot);
  const int dimv = robot.dimv();
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  const int dimx = 2 * robot.dimv();
  seed = Eigen::MatrixXd::Random(dimx, dimx);
  riccati.N = seed * seed.transpose();
  riccati.Pi.setRandom();
  riccati.pi.setRandom();
  return riccati; 
}


void ImpulseRiccatiFactorizerTest::testBackwardRecursion(const Robot& robot) {
  const int dimv = robot.dimv();
  const RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  ImpulseKKTMatrix kkt_matrix = createKKTMatrix(robot);
  ImpulseKKTResidual kkt_residual = createKKTResidual(robot);
  ImpulseKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseRiccatiFactorizer factorizer(robot);
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  RiccatiFactorization riccati(robot), riccati_ref(robot);
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati_ref.Pqq));
  EXPECT_TRUE(riccati.Pqv.isApprox(riccati_ref.Pqv));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati_ref.Pvq));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati_ref.Pvv));
  EXPECT_TRUE(riccati.sq.isApprox(riccati_ref.sq));
  EXPECT_TRUE(riccati.sv.isApprox(riccati_ref.sv));
  EXPECT_TRUE(riccati.Pqq.isApprox(riccati.Pqq.transpose()));
  EXPECT_TRUE(riccati.Pvv.isApprox(riccati.Pvv.transpose()));
  EXPECT_TRUE(riccati.Pvq.isApprox(riccati.Pqv.transpose()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
}


void ImpulseRiccatiFactorizerTest::testForwardRecursion(const Robot& robot) {
  const int dimv = robot.dimv();
  RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  RiccatiFactorization riccati_next_ref = riccati_next;
  ImpulseKKTMatrix kkt_matrix = createKKTMatrix(robot);
  ImpulseKKTResidual kkt_residual = createKKTResidual(robot);
  ImpulseKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseRiccatiFactorizer factorizer(robot);
  RiccatiFactorization riccati(robot);
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  factorizer.forwardRiccatiRecursionSerial(riccati, kkt_matrix, kkt_residual, riccati_next);
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_ref(robot);
  RiccatiFactorization riccati_ref(robot);
  backward_recursion_ref.factorizeKKTMatrix(riccati_next_ref, kkt_matrix_ref, kkt_residual_ref);
  backward_recursion_ref.factorizeRiccatiFactorization(riccati_next_ref, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
  if (!robot.has_floating_base()) kkt_matrix_ref.Fqq().setIdentity();
  riccati_next_ref.Pi = kkt_matrix_ref.Fxx() * riccati_ref.Pi;
  riccati_next_ref.pi = kkt_residual_ref.Fx() + kkt_matrix_ref.Fxx() * riccati_ref.pi;
  riccati_next_ref.N = kkt_matrix_ref.Fxx() * riccati_ref.N * kkt_matrix_ref.Fxx().transpose();
  EXPECT_TRUE(riccati_next.Pi.isApprox(riccati_next_ref.Pi));
  EXPECT_TRUE(riccati_next.pi.isApprox(riccati_next_ref.pi));
  EXPECT_TRUE(riccati_next.N.isApprox(riccati_next_ref.N));
  riccati.n.setRandom();
  const Eigen::VectorXd dx0 = Eigen::VectorXd::Random(2*dimv);
  ImpulseSplitDirection d(robot), d_ref(robot);
  factorizer.computeStateDirection(riccati, dx0, d);
  d_ref.dx().noalias() = riccati.Pi * dx0 + riccati.pi - riccati.N * riccati.n;
  EXPECT_TRUE(d.isApprox(d_ref));
  factorizer.computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq + riccati.n.head(dimv);
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv + riccati.n.tail(dimv);
  EXPECT_TRUE(d.isApprox(d_ref));
}


void ImpulseRiccatiFactorizerTest::testFactorizeStateConstraintFactorization(const Robot& robot) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  const RiccatiFactorization riccati_next = createRiccatiFactorization(robot);
  const ImpulseKKTMatrix kkt_matrix = createKKTMatrix(robot);
  const ImpulseKKTResidual kkt_residual = createKKTResidual(robot);
  ImpulseRiccatiFactorizer factorizer(robot);
  RiccatiFactorization riccati(robot), riccati_ref(robot);
  const int dimf = 6;
  Eigen::MatrixXd T_next(2*dimv, dimf), T(2*dimv, dimf), T_ref(2*dimv, dimf);
  T_next.setRandom();
  T.setZero();
  T_ref.setZero();
  factorizer.backwardStateConstraintFactorization(T_next, kkt_matrix, T);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  if (robot.has_floating_base()) {
    A.topLeftCorner(dimv, dimv) = kkt_matrix.Fqq();
  }
  else {
    A.topLeftCorner(dimv, dimv).setIdentity();
  }
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix.Fvv();
  T_ref = A.transpose() * T_next;
  EXPECT_TRUE(T.isApprox(T_ref));
}


TEST_F(ImpulseRiccatiFactorizerTest, fixedBase) {
  testBackwardRecursion(fixed_base_robot);
  testForwardRecursion(fixed_base_robot);
  testFactorizeStateConstraintFactorization(fixed_base_robot);
}


TEST_F(ImpulseRiccatiFactorizerTest, floating_base) {
  testBackwardRecursion(floating_base_robot);
  testForwardRecursion(floating_base_robot);
  testFactorizeStateConstraintFactorization(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}