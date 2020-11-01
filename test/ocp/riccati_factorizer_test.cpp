#include <string>
#include <memory>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/ocp/riccati_gain.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"


namespace idocp {

class RiccatiFactorizerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    fixed_base_urdf = "../urdf/iiwa14/iiwa14.urdf";
    floating_base_urdf = "../urdf/anymal/anymal.urdf";
    fixed_base_robot = Robot(fixed_base_urdf);
    floating_base_robot = Robot(floating_base_urdf);
    dtau = std::abs(Eigen::VectorXd::Random(1)[0]);
  }

  virtual void TearDown() {
  }

  static void testBackwardRecursionUnits(const Robot& robot, const double dtau);

  static void testBackwardRecursion(const Robot& robot, const double dtau);

  static void testForwardRecursion(const Robot& robot, const double dtau);

  static void testComputeCostateDirection(const Robot& robot, const double dtau);

  static void testComputeControlInputDirection(const Robot& robot, const double dtau);

  double dtau;
  std::string fixed_base_urdf, floating_base_urdf;
  Robot fixed_base_robot, floating_base_robot;
};


void RiccatiFactorizerTest::testBackwardRecursionUnits(const Robot& robot, const double dtau) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiSolution riccati_next(robot);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
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
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd Qqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Qvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Quu = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Qqu() = Qqu;
  kkt_matrix.Qvu() = Qvu;
  kkt_matrix.Quu() = Quu;
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  const KKTMatrix kkt_matrix_ref = kkt_matrix;
  const KKTResidual kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot);
  factorizer.factorizeMatrices(riccati_next, dtau, kkt_matrix, kkt_residual);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  A.topLeftCorner(dimv, dimv) = kkt_matrix_ref.Fqq();
  A.topRightCorner(dimv, dimv) = kkt_matrix_ref.Fqv();
  A.bottomLeftCorner(dimv, dimv) = kkt_matrix_ref.Fvq();
  A.bottomRightCorner(dimv, dimv) = kkt_matrix_ref.Fvv();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*dimv, dimu);
  B.bottomRows(dimv) = kkt_matrix_ref.Fvu();
  Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(2*dimv, 2*dimv);
  P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  const Eigen::MatrixXd F_ref = kkt_matrix_ref.Qxx() + A.transpose() * P_next * A;
  const Eigen::MatrixXd H_ref = kkt_matrix_ref.Qxu() + A.transpose() * P_next * B;
  const Eigen::MatrixXd G_ref = kkt_matrix_ref.Quu() + B.transpose() * P_next * B;
  Eigen::VectorXd sx_next = Eigen::VectorXd::Zero(2*dimv);
  sx_next.head(dimv) = riccati_next.sq;
  sx_next.tail(dimv) = riccati_next.sv;
  const Eigen::VectorXd lu_ref = B.transpose() * P_next * kkt_residual_ref.Fx() - B.transpose() * sx_next + kkt_residual_ref.lu();
  EXPECT_TRUE(F_ref.isApprox(kkt_matrix.Qxx()));
  EXPECT_TRUE(kkt_matrix.Qxx().isApprox(kkt_matrix.Qxx().transpose()));
  EXPECT_TRUE(H_ref.isApprox(kkt_matrix.Qxu()));
  EXPECT_TRUE(G_ref.isApprox(kkt_matrix.Quu()));
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
  EXPECT_TRUE(lu_ref.isApprox(kkt_residual.lu()));
  RiccatiGain gain(robot), gain_ref(robot);
  factorizer.computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual, gain);
  const Eigen::MatrixXd Ginv = G_ref.inverse();
  gain_ref.K = - Ginv * H_ref.transpose();
  gain_ref.k = - Ginv * lu_ref;
  EXPECT_TRUE(gain.K.isApprox(gain_ref.K));
  EXPECT_TRUE(gain.k.isApprox(gain_ref.k));
  RiccatiSolution riccati(robot);
  factorizer.factorizeRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, gain, riccati);
  // const Eigen::MatrixXd P_ref = F_ref + H_ref * gain.K;
  const Eigen::MatrixXd P_ref = F_ref - gain.K.transpose() * G_ref * gain.K;
  const Eigen::VectorXd s_ref = A.transpose() * sx_next - A.transpose() * P_next * kkt_residual_ref.Fx() - kkt_residual_ref.lx() - H_ref * gain.k;
  EXPECT_TRUE(P_ref.topLeftCorner(dimv, dimv).isApprox(riccati.Pqq));
  EXPECT_TRUE(P_ref.topRightCorner(dimv, dimv).isApprox(riccati.Pqv));
  EXPECT_TRUE(P_ref.bottomLeftCorner(dimv, dimv).isApprox(riccati.Pvq));
  EXPECT_TRUE(P_ref.bottomRightCorner(dimv, dimv).isApprox(riccati.Pvv));
  EXPECT_TRUE(s_ref.head(dimv).isApprox(riccati.sq));
  EXPECT_TRUE(s_ref.tail(dimv).isApprox(riccati.sv));
}


void RiccatiFactorizerTest::testBackwardRecursion(const Robot& robot, const double dtau) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiSolution riccati_next(robot);
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
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
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd Qqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Qvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Quu = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Qqu() = Qqu;
  kkt_matrix.Qvu() = Qvu;
  kkt_matrix.Quu() = Quu;
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  KKTMatrix kkt_matrix_ref = kkt_matrix;
  KKTResidual kkt_residual_ref = kkt_residual;
  RiccatiFactorizer factorizer(robot), factorizer_ref(robot);
  RiccatiSolution riccati(robot), riccati_ref(robot);
  RiccatiGain gain(robot), gain_ref(robot);
  factorizer.factorizeBackwardRicursion(riccati_next, dtau, kkt_matrix, kkt_residual, gain, riccati);
  factorizer_ref.factorizeMatrices(riccati_next, dtau, kkt_matrix_ref, kkt_residual_ref);
  factorizer_ref.computeFeedbackGainAndFeedforward(kkt_matrix_ref, kkt_residual_ref, gain_ref);
  factorizer_ref.factorizeRecursion(riccati_next, dtau, kkt_matrix_ref, kkt_residual_ref, gain_ref, riccati_ref);
  EXPECT_TRUE(gain.K.isApprox(gain_ref.K));
  EXPECT_TRUE(gain.k.isApprox(gain_ref.k));
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
  EXPECT_TRUE(kkt_matrix.Quu().isApprox(kkt_matrix.Quu().transpose()));
  std::cout << riccati.Pqq - riccati.Pqq.transpose() << std::endl;
  std::cout << riccati.Pqv - riccati.Pvq.transpose() << std::endl;
  std::cout << riccati.Pvv - riccati.Pvv.transpose() << std::endl;
  std::cout << kkt_matrix.Qxx() - kkt_matrix.Qxx().transpose() << std::endl;
  std::cout << kkt_matrix.Quu() - kkt_matrix.Quu().transpose() << std::endl;
}


void RiccatiFactorizerTest::testForwardRecursion(const Robot& robot, const double dtau) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  KKTMatrix kkt_matrix(robot);
  KKTResidual kkt_residual(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qqq = seed * seed.transpose();
  const Eigen::MatrixXd Qqv = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvq = Qqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  const Eigen::MatrixXd Qvv = seed * seed.transpose();
  seed = Eigen::MatrixXd::Random(dimu, dimu);
  const Eigen::MatrixXd Qqu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Qvu = Eigen::MatrixXd::Random(dimv, dimu);
  const Eigen::MatrixXd Quu = seed * seed.transpose();
  kkt_matrix.Qqq() = Qqq;
  kkt_matrix.Qqv() = Qqv;
  kkt_matrix.Qvq() = Qvq;
  kkt_matrix.Qvv() = Qvv;
  kkt_matrix.Qqu() = Qqu;
  kkt_matrix.Qvu() = Qvu;
  kkt_matrix.Quu() = Quu;
  kkt_matrix.Fqq() = Eigen::MatrixXd::Identity(dimv, dimv);
  if (robot.has_floating_base()) {
    kkt_matrix.Fqq().topLeftCorner(6, 6).setRandom();
  }
  kkt_matrix.Fqv() = dtau * Eigen::MatrixXd::Identity(dimv, dimv);
  kkt_matrix.Fvq().setRandom();
  kkt_matrix.Fvv().setRandom();
  kkt_matrix.Fvu().setRandom();
  kkt_residual.lu().setRandom();
  kkt_residual.Fx().setRandom();
  const SplitDirection d = SplitDirection::Random(robot);
  SplitDirection d_next(robot), d_next_ref(robot);
  RiccatiFactorizer factorizer(robot);
  factorizer.factorizeForwardRicursion(kkt_matrix, kkt_residual, d, dtau, d_next);
  d_next_ref.dx() = kkt_matrix.Fxx() * d.dx() + kkt_matrix.Fxu() * d.du() + kkt_residual.Fx();
  EXPECT_TRUE(d_next.isApprox(d_next_ref));
}


void RiccatiFactorizerTest::testComputeCostateDirection(const Robot& robot, const double dtau) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiSolution riccati(robot);
  Eigen::MatrixXd seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pqq = seed * seed.transpose();
  riccati.Pqv = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvq = riccati.Pqv.transpose();
  seed = Eigen::MatrixXd::Random(dimv, dimv);
  riccati.Pvv = seed * seed.transpose();
  riccati.sq.setRandom();
  riccati.sv.setRandom();
  SplitDirection d = SplitDirection::Random(robot);
  d.dlmd().setZero();
  d.dgmm().setZero();
  SplitDirection d_ref = d;
  RiccatiFactorizer::computeCostateDirection(riccati, d);
  d_ref.dlmd() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d_ref.dgmm() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
  EXPECT_TRUE(d.isApprox(d_ref));
}


void RiccatiFactorizerTest::testComputeControlInputDirection(const Robot& robot, const double dtau) {
  const int dimv = robot.dimv();
  const int dimu = robot.dimu();
  RiccatiGain gain(robot);
  gain.K.setRandom();
  gain.k.setRandom();
  SplitDirection d = SplitDirection::Random(robot);
  d.du().setZero();
  SplitDirection d_ref = d;
  RiccatiFactorizer::computeControlInputDirection(gain, d);
  d_ref.du() = gain.K * d_ref.dx() + gain.k;
  EXPECT_TRUE(d.isApprox(d_ref));
}


TEST_F(RiccatiFactorizerTest, fixedBase) {
  testBackwardRecursionUnits(fixed_base_robot, dtau);
  testBackwardRecursion(fixed_base_robot, dtau);
  testForwardRecursion(fixed_base_robot, dtau);
  testComputeCostateDirection(fixed_base_robot, dtau);
  testComputeControlInputDirection(fixed_base_robot, dtau);
}


TEST_F(RiccatiFactorizerTest, floating_base) {
  testBackwardRecursionUnits(floating_base_robot, dtau);
  testBackwardRecursion(floating_base_robot, dtau);
  testForwardRecursion(floating_base_robot, dtau);
  testComputeCostateDirection(floating_base_robot, dtau);
  testComputeControlInputDirection(floating_base_robot, dtau);
}

} // namespace idocp


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}