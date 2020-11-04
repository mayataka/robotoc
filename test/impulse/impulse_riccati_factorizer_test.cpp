#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_solution.hpp"
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

  static void testBackwardRecursionUnits(const Robot& robot);
  static void testBackwardRecursion(const Robot& robot);
  static void testForwardRecursion(const Robot& robot);
  static void testComputeCostateDirection(const Robot& robot);
  static void testComputeControlInputDirection(const Robot& robot);

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


void ImpulseRiccatiFactorizerTest::testBackwardRecursionUnits(const Robot& robot) {
  const int dimv = robot.dimv();
  RiccatiSolution riccati_next(robot);
  ImpulseKKTMatrix kkt_matrix(robot);
  ImpulseKKTResidual kkt_residual(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pqq = seed * seed.transpose();
  riccati_next.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvv = seed * seed.transpose();
  riccati_next.sq.setRandom();
  riccati_next.sv.setRandom();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = seed * seed.transpose();
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Qqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Fqq().setIdentity();
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_residual.Fx().setRandom();
  const ImpulseKKTMatrix kkt_matrix_ref = kkt_matrix;
  const ImpulseKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseRiccatiFactorizer factorizer(robot);
  factorizer.factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Fqq();
  A.topRightCorner(dimv, dimv).setZero();
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix_ref.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix_ref.Fvv();
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx() + A.transpose() * P_next * A;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  RiccatiSolution riccati(robot);
  factorizer.factorizeRiccatiSolution(riccati_next, kkt_matrix, kkt_residual, riccati);
  const Eigen::MatrixXd P_ref = F_ref;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * kkt_residual_ref.Fx() - kkt_residual_ref.lx();
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}


void ImpulseRiccatiFactorizerTest::testBackwardRecursion(const Robot& robot) {
  const int dimv = robot.dimv();
  RiccatiSolution riccati_next(robot);
  ImpulseKKTMatrix kkt_matrix(robot);
  ImpulseKKTResidual kkt_residual(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pqq = seed * seed.transpose();
  riccati_next.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvq = riccati_next.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati_next.Pvv = seed * seed.transpose();
  riccati_next.sq.setRandom();
  riccati_next.sv.setRandom();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = seed * seed.transpose();
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Qqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Fqq().setIdentity(); 
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_residual.Fx().setRandom();
  ImpulseKKTMatrix kkt_matrix_ref = kkt_matrix;
  ImpulseKKTResidual kkt_residual_ref = kkt_residual;
  ImpulseRiccatiFactorizer factorizer(robot), factorizer_ref(robot);
  RiccatiSolution riccati(robot), riccati_ref(robot);
  factorizer.backwardRiccatiRecursion(riccati_next, kkt_matrix, kkt_residual, riccati);
  factorizer_ref.factorizeKKTMatrix(riccati_next, kkt_matrix_ref, kkt_residual_ref);
  factorizer_ref.factorizeRiccatiSolution(riccati_next, kkt_matrix_ref, kkt_residual_ref, riccati_ref);
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
  ImpulseKKTMatrix kkt_matrix(robot);
  ImpulseKKTResidual kkt_residual(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = seed * seed.transpose();
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Qqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.has_floating_base()) {
    kkt_matrix.Fqq().topLeftCorner(6, 6).setRandom();
  }
  kkt_matrix.Fqv().setZero();
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_residual.Fx().setRandom();
  const ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot);
  SplitDirection d_next(robot), d_next_ref(robot);
  ImpulseRiccatiFactorizer factorizer(robot);
  factorizer.forwardRiccatiRecursion(kkt_matrix, kkt_residual, d, d_next);
  d_next_ref.dx() = kkt_matrix.Fxx() * d.dx() + kkt_residual.Fx();
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
}


void ImpulseRiccatiFactorizerTest::testComputeCostateDirection(const Robot& robot) {
  const int dimv = robot.dimv();
  RiccatiSolution riccati(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  ImpulseSplitDirection d = ImpulseSplitDirection::Random(robot);
  d.dlmd().setZero();
  d.dgmm().setZero();
  ImpulseSplitDirection d_ref = d;
  ImpulseRiccatiFactorizer::computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(ImpulseRiccatiFactorizerTest, fixedBase) {
  testBackwardRecursionUnits(fixed_base_robot);
  testBackwardRecursion(fixed_base_robot);
  testForwardRecursion(fixed_base_robot);
  testComputeCostateDirection(fixed_base_robot);
}


TEST_F(ImpulseRiccatiFactorizerTest, floating_base) {
  testBackwardRecursionUnits(floating_base_robot);
  testBackwardRecursion(floating_base_robot);
  testForwardRecursion(floating_base_robot);
  testComputeCostateDirection(floating_base_robot);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}